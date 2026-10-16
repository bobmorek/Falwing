/* FAL Wing swell field: linear-wave refraction, shoaling and depth-limited breaking over real bathymetry.
 *
 * 1. Grid (~75 m) over the coastline box. Land from the OSM coastline, depth from the EMODnet DTM.
 * 2. Dispersion: k from w^2 = g k tanh(kh) (Fenton & McKee explicit form), phase speed c = w/k,
 *    group speed cg = c/2 (1 + 2kh / sinh 2kh).
 * 3. Wave-front arrival time tau from the eikonal equation |grad tau| = 1/c, solved by fast marching
 *    from the open-sea boundary with a plane wave. Crests are lines of equal tau (mod period), so they
 *    slow and bend toward shallow water (refraction) and wrap round headlands (diffraction-like).
 * 4. Energy along rays: wave-ray spacing b grows as d(ln b)/ds = div(u), u = unit propagation
 *    direction, so the refraction coefficient Kr = sqrt(b0/b) = exp(-1/2 * integral div(u) ds).
 *    Shoaling Ks = sqrt(cg0/cg). Local height H = Hs*Kr*Ks, capped at 0.78*h (depth-limited breaking).
 * Illustrative only: no wind-sea growth, currents, tides, bottom friction or proper diffraction.
 */
"use strict";
(function(){
  const g = 9.81, CELL = 75;
  let G = null;            // grid
  let F = null;            // current solved field
  let SW = null;           // current swell {hs,tp,dir}


  function kh_of(T, h){                       // Fenton & McKee (1990)
    const w = 2*Math.PI/T, kh0 = w*w/g*h;
    return kh0 / Math.pow(Math.tanh(Math.pow(kh0, 0.75)), 2/3);
  }
  function bilinear(a, nx, ny, fx, fy){
    fx = Math.max(0, Math.min(nx-1.001, fx)); fy = Math.max(0, Math.min(ny-1.001, fy));
    const i = fx|0, j = fy|0, dx = fx-i, dy = fy-j, k = j*nx+i;
    return a[k]*(1-dx)*(1-dy) + a[k+1]*dx*(1-dy) + a[k+nx]*(1-dx)*dy + a[k+nx+1]*dx*dy;
  }

  function build(coast, bathy){
    const [W, S, E, N] = coast.bbox, lat0 = (S+N)/2;
    const dlon = CELL/(111320*Math.cos(lat0*Math.PI/180)), dlat = CELL/110574;
    const nx = Math.ceil((E-W)/dlon), ny = Math.ceil((N-S)/dlat), n = nx*ny;

    // land mask from the coastline polygons
    const cv = document.createElement("canvas"); cv.width = nx; cv.height = ny;
    const cx = cv.getContext("2d"); cx.fillStyle = "#000"; cx.beginPath();
    const polys = coast.land.type === "MultiPolygon" ? coast.land.coordinates : [coast.land.coordinates];
    for (const poly of polys) for (const ring of poly) {
      ring.forEach(([lo, la], k) => { const x = (lo-W)/dlon, y = (N-la)/dlat; k ? cx.lineTo(x, y) : cx.moveTo(x, y); });
      cx.closePath();
    }
    cx.fill("evenodd");
    const px = cx.getImageData(0, 0, nx, ny).data, land = new Uint8Array(n);
    for (let k = 0; k < n; k++) land[k] = px[k*4+3] > 110 ? 1 : 0;

    // distance to land (chamfer), metres
    const dist = new Float32Array(n), BIG = 1e9, R2 = Math.SQRT2;
    for (let k = 0; k < n; k++) dist[k] = land[k] ? 0 : BIG;
    for (let j = 0; j < ny; j++) for (let i = 0; i < nx; i++) { const k = j*nx+i; let d = dist[k]; if (!d) continue;
      if (i > 0) d = Math.min(d, dist[k-1]+1);
      if (j > 0) { d = Math.min(d, dist[k-nx]+1); if (i > 0) d = Math.min(d, dist[k-nx-1]+R2); if (i < nx-1) d = Math.min(d, dist[k-nx+1]+R2); }
      dist[k] = d; }
    for (let j = ny-1; j >= 0; j--) for (let i = nx-1; i >= 0; i--) { const k = j*nx+i; let d = dist[k]; if (!d) continue;
      if (i < nx-1) d = Math.min(d, dist[k+1]+1);
      if (j < ny-1) { d = Math.min(d, dist[k+nx]+1); if (i < nx-1) d = Math.min(d, dist[k+nx+1]+R2); if (i > 0) d = Math.min(d, dist[k+nx-1]+R2); }
      dist[k] = d; }
    for (let k = 0; k < n; k++) dist[k] *= CELL;

    // depth: EMODnet where available, otherwise a gentle slope from the shore
    let bd = null;
    if (bathy && bathy.depth_b64) { const s = atob(bathy.depth_b64); bd = new Uint8Array(s.length); for (let i = 0; i < s.length; i++) bd[i] = s.charCodeAt(i); }
    const depth = new Float32Array(n);
    for (let j = 0; j < ny; j++) for (let i = 0; i < nx; i++) {
      const k = j*nx+i; if (land[k]) continue;
      const lo = W+(i+.5)*dlon, la = N-(j+.5)*dlat, slope = Math.min(45, 1+0.025*dist[k]);
      let h = NaN;
      if (bd) { const [bW, bS, bE, bN] = bathy.bbox;
        h = bilinear(bd, bathy.nx, bathy.ny, (lo-bW)/(bE-bW)*bathy.nx-.5, (bN-la)/(bN-bS)*bathy.ny-.5)*bathy.scale; }
      depth[k] = h > 0.5 ? h : (bd ? Math.min(slope, 2) : slope);   // DTM says dry but coastline says water: very shallow
    }
    const coastal = new Uint8Array(n);
    for (let k = 0; k < n; k++) coastal[k] = dist[k] <= 1.6*CELL ? 1 : 0;   // includes the extended land rim
    return { W, N, nx, ny, n, dlon, dlat, land, depth, coastal, hasBathy: !!bd, cache: new Map() };
  }

  function solve(T, dir){
    T = Math.max(3, Math.round(T*2)/2); dir = ((Math.round(dir/5)*5)%360+360)%360;
    const key = dir+"|"+T; if (G.cache.has(key)) return G.cache.get(key);
    const { nx, ny, n, land, depth } = G, H = CELL, w = 2*Math.PI/T, c0 = g*T/(2*Math.PI), cg0 = c0/2;

    const slow = new Float32Array(n), Ks = new Float32Array(n); const slowAvg = 1/c0;
    for (let k = 0; k < n; k++) { if (land[k]) continue;
      const h = depth[k], kh = kh_of(T, h), c = w*h/kh, cg = kh > 15 ? c/2 : c/2*(1+2*kh/Math.sinh(2*kh));
      slow[k] = 1/c; Ks[k] = Math.sqrt(cg0/cg); }

    // propagation direction in grid axes (i east, j south)
    const b = (dir+180)*Math.PI/180, vi = Math.sin(b), vj = -Math.cos(b);
    const pmin = Math.min(0, (nx-1)*vi, (ny-1)*vj, (nx-1)*vi+(ny-1)*vj);

    const tau = new Float32Array(n).fill(Infinity), state = new Uint8Array(n), order = new Int32Array(n), src = new Uint8Array(n);
    let on = 0;
    // binary min-heap with lazy deletion
    let heap = new Int32Array(n*5), hn = 0;
    const push = k => { let i = hn++; heap[i] = k; while (i) { const p = (i-1)>>1; if (tau[heap[p]] <= tau[heap[i]]) break; const t = heap[p]; heap[p] = heap[i]; heap[i] = t; i = p; } };
    const pop = () => { const top = heap[0]; heap[0] = heap[--hn]; let i = 0;
      for (;;) { const l = 2*i+1, r = l+1; let m = i;
        if (l < hn && tau[heap[l]] < tau[heap[m]]) m = l; if (r < hn && tau[heap[r]] < tau[heap[m]]) m = r;
        if (m === i) break; const t = heap[m]; heap[m] = heap[i]; heap[i] = t; i = m; }
      return top; };

    const seed = (i, j, ni, nj) => { const k = j*nx+i; if (land[k] || ni*vi+nj*vj <= 0.05 || src[k]) return;
      tau[k] = (i*vi+j*vj-pmin)*H/c0; state[k] = 1; src[k] = 1; push(k); };
    for (let i = 0; i < nx; i++) { seed(i, 0, 0, 1); seed(i, ny-1, 0, -1); }
    for (let j = 0; j < ny; j++) { seed(0, j, 1, 0); seed(nx-1, j, -1, 0); }

    const known = m => state[m] === 2 ? tau[m] : Infinity;
    const update = m => { const i = m % nx, j = (m/nx)|0;
      let a = Math.min(i > 0 ? known(m-1) : Infinity, i < nx-1 ? known(m+1) : Infinity);
      let c = Math.min(j > 0 ? known(m-nx) : Infinity, j < ny-1 ? known(m+nx) : Infinity);
      if (a > c) { const t = a; a = c; c = t; }
      const f = slow[m]*H;
      if (c === Infinity || c-a >= f) return a+f;
      return (a+c+Math.sqrt(2*f*f-(c-a)*(c-a)))/2; };

    while (hn) { const k = pop(); if (state[k] === 2) continue; state[k] = 2; order[on++] = k;
      const i = k % nx, j = (k/nx)|0;
      const nb = [i > 0 ? k-1 : -1, i < nx-1 ? k+1 : -1, j > 0 ? k-nx : -1, j < ny-1 ? k+nx : -1];
      for (const m of nb) { if (m < 0 || land[m] || state[m] === 2) continue;
        const t = update(m); if (t < tau[m]) { tau[m] = t; state[m] = 1; if (hn >= heap.length-1) { const h2 = new Int32Array(heap.length*2); h2.set(heap); heap = h2; } push(m); } } }

    // unit propagation direction u = grad(tau)/|grad(tau)|
    const ux = new Float32Array(n), uy = new Float32Array(n);
    const fin = m => m >= 0 && m < n && !land[m] && tau[m] < Infinity;
    for (let k = 0; k < n; k++) { if (!fin(k)) continue; const i = k % nx, j = (k/nx)|0;
      const l = i > 0 && fin(k-1), r = i < nx-1 && fin(k+1), u = j > 0 && fin(k-nx), d = j < ny-1 && fin(k+nx);
      const gx = l && r ? (tau[k+1]-tau[k-1])/2 : r ? tau[k+1]-tau[k] : l ? tau[k]-tau[k-1] : 0;
      const gy = u && d ? (tau[k+nx]-tau[k-nx])/2 : d ? tau[k+nx]-tau[k] : u ? tau[k]-tau[k-nx] : 0;
      const m = Math.hypot(gx, gy); if (m > 0) { ux[k] = gx/m; uy[k] = gy/m; } else { ux[k] = vi; uy[k] = vj; } }

    // ray spreading integrated along the march order: L = integral div(u) ds
    const L = new Float32Array(n);
    for (let o = 0; o < on; o++) { const k = order[o]; if (src[k]) { L[k] = 0; continue; }
      const i = k % nx, j = (k/nx)|0;
      let div = 0;
      div += (i > 0 && i < nx-1 && fin(k-1) && fin(k+1)) ? (ux[k+1]-ux[k-1])/2 : 0;
      div += (j > 0 && j < ny-1 && fin(k-nx) && fin(k+nx)) ? (uy[k+nx]-uy[k-nx])/2 : 0;
      div = Math.max(-0.35, Math.min(0.35, div));
      let sw = 0, sl = 0;
      const up = (m, di, dj) => { if (!fin(m) || tau[m] >= tau[k]) return; const wgt = Math.max(0, -(di*ux[k]+dj*uy[k])); sw += wgt; sl += wgt*L[m]; };
      if (i > 0) up(k-1, -1, 0); if (i < nx-1) up(k+1, 1, 0); if (j > 0) up(k-nx, 0, -1); if (j < ny-1) up(k+nx, 0, 1);
      L[k] = (sw > 0 ? sl/sw : 0) + div;
    }

    const rel = new Float32Array(n);
    for (let k = 0; k < n; k++) { if (!fin(k)) continue;
      const Kr = Math.max(0.03, Math.min(2.2, Math.exp(-L[k]/2)));
      const turn = ux[k]*vi + uy[k]*vj;                      // cos of total turning
      const sh0 = Math.min(1, Math.max(0, (1+turn)/1.4)), shadow = sh0*sh0;   // little energy survives turning past ~100 deg
      rel[k] = Kr*Math.min(Ks[k], 2.5)*shadow; }

    // extend the field 2 cells under the land so drawn crests run right up to the coastline
    const tx = Float32Array.from(tau), rx = Float32Array.from(rel);
    for (let pass = 0; pass < 2; pass++) {
      const t0 = Float32Array.from(tx), r0 = Float32Array.from(rx);
      for (let k = 0; k < n; k++) { if (t0[k] < Infinity) continue; const i = k % nx, j = (k/nx)|0;
        let st = 0, sr = 0, c = 0;
        const add = m => { if (t0[m] < Infinity) { st += t0[m]; sr += r0[m]; c++; } };
        if (i > 0) add(k-1); if (i < nx-1) add(k+1); if (j > 0) add(k-nx); if (j < ny-1) add(k+nx);
        if (c) { tx[k] = st/c + 0.15*slowAvg*H; rx[k] = sr/c; } } }
    // light smoothing of tau removes the first-order marching staircase without moving the fronts
    for (let pass = 0; pass < 2; pass++) { const t0 = Float32Array.from(tx);
      for (let j = 1; j < ny-1; j++) for (let i = 1; i < nx-1; i++) { const k = j*nx+i; if (!(t0[k] < Infinity)) continue;
        let s = 0, c = 0;
        for (let dj = -1; dj <= 1; dj++) for (let di = -1; di <= 1; di++) { const v = t0[k+dj*nx+di]; if (v < Infinity) { const w = (di||dj) ? 1 : 2; s += w*v; c += w; } }
        tx[k] = s/c; } }
    const out = { T, dir, tau: tx, rel: rx };
    if (G.cache.size > 24) G.cache.delete(G.cache.keys().next().value);
    G.cache.set(key, out);
    return out;
  }

  // per-cell crest and foam strength for the current height
  function shade(){
    if (!G || !F || !SW) return;
    const { n, depth, coastal } = G, Hs = Math.max(0.05, SW.hs);
    const crest = F.crest || (F.crest = new Float32Array(n)), foam = F.foam || (F.foam = new Float32Array(n));
    for (let k = 0; k < n; k++) { const r = F.rel[k]; if (!r) { crest[k] = 0; foam[k] = 0; continue; }
      const Hl = Math.min(Hs*r, 0.78*Math.max(depth[k], 0.6));
      const rr = Math.min(1, Hl/Math.max(0.3, Hs)); crest[k] = rr < 0.12 ? 0 : Math.pow(rr, 0.8) * (0.6 + 0.4*Math.min(Hs/2, 1));
      foam[k] = coastal[k] && Hl > 0.15 ? Math.min(1, Hl/0.6) : 0; }
  }

  let pending = 0;
  function resolve(){
    if (!G || !SW || SW.dir == null || !SW.tp) { F = null; return; }
    clearTimeout(pending);
    pending = setTimeout(() => { F = solve(SW.tp, SW.dir); shade(); }, 30);
  }

  // Draw crests as isolines of phase = (tau - t)/Tshow = integer, by marching squares over the visible
  // cells. tau is continuous, so each integer level is one crest. Strength sets opacity; near-shore
  // crests are drawn as foam.
  const NB = 6;
  function draw(ctx, map, phase, Tshow, col, foamCol){
    if (!G || !F || !F.crest) return false;
    const { nx, ny, W, N, dlon, dlat } = G, tau = F.tau, crest = F.crest, foam = F.foam;
    const b = map.getBounds();
    const i0 = Math.max(0, Math.floor((b.getWest()-W)/dlon)-1), i1 = Math.min(nx-1, Math.ceil((b.getEast()-W)/dlon)+1);
    const j0 = Math.max(0, Math.floor((N-b.getNorth())/dlat)-1), j1 = Math.min(ny-1, Math.ceil((N-b.getSouth())/dlat)+1);
    if (i1 <= i0 || j1 <= j0) return true;
    // screen coords of cell centres, unrounded (Leaflet's latLngToContainerPoint rounds to whole pixels):
    // x is linear in lon, y is computed per row (Mercator)
    const zm = map.getZoom(), po = map.getPixelOrigin(), o = map.containerPointToLayerPoint([0, 0]);
    const P = (la, lo) => { const p = map.project([la, lo], zm); return [p.x-po.x-o.x, p.y-po.y-o.y]; };
    const pa = { x: P(N, W+0.5*dlon)[0] }, sx = P(N, W+1.5*dlon)[0]-pa.x;
    const ys = new Float64Array(ny); for (let j = j0; j <= j1; j++) ys[j] = P(N-(j+.5)*dlat, W)[1];
    const X = i => pa.x + i*sx;
    const inv = 1/Tshow, ph = phase*inv;
    const cp = [], fp = []; for (let q = 0; q < NB; q++) { cp.push(new Path2D()); fp.push(new Path2D()); }
    const v = [0,0,0,0], ex = [0,0,0,0], ey = [0,0,0,0];
    for (let j = j0; j < j1; j++) { const y0 = ys[j], y1 = ys[j+1];
      for (let i = i0; i < i1; i++) {
        const k = j*nx+i, a = tau[k], bb = tau[k+1], c = tau[k+nx+1], d = tau[k+nx];
        if (!(a < Infinity && bb < Infinity && c < Infinity && d < Infinity)) continue;
        const st = (crest[k]+crest[k+1]+crest[k+nx]+crest[k+nx+1])*0.25;
        const fo = Math.max(foam[k], foam[k+1], foam[k+nx], foam[k+nx+1]);
        if (st < 0.06 && fo < 0.06) continue;
        v[0] = a*inv-ph; v[1] = bb*inv-ph; v[2] = c*inv-ph; v[3] = d*inv-ph;   // corners: TL, TR, BR, BL
        const lo = Math.min(v[0], v[1], v[2], v[3]), hi = Math.max(v[0], v[1], v[2], v[3]);
        let lev = Math.ceil(lo); if (lev > hi) continue;
        const xL = X(i), xR = X(i+1);
        const cx = [xL, xR, xR, xL], cy = [y0, y0, y1, y1];
        for (; lev <= hi; lev++) {
          let m = 0;
          for (let e = 0; e < 4; e++) { const p = v[e]-lev, q = v[(e+1)&3]-lev;
            if ((p < 0) !== (q < 0)) { const t = p/(p-q), e2 = (e+1)&3;
              ex[m] = cx[e]+(cx[e2]-cx[e])*t; ey[m] = cy[e]+(cy[e2]-cy[e])*t; m++; } }
          if (m < 2) continue;
          const useFoam = fo > 0.08 && fo*1.1 > st;
          const path = useFoam ? fp[Math.min(NB-1, (fo*NB)|0)] : cp[Math.min(NB-1, (st*NB)|0)];
          path.moveTo(ex[0], ey[0]); path.lineTo(ex[1], ey[1]);
          if (m === 4) { path.moveTo(ex[2], ey[2]); path.lineTo(ex[3], ey[3]); }
        }
      } }
    const z = map.getZoom(), lw = Math.max(1.35, Math.min(2.6, 0.6 + (z-10)*0.4));
    ctx.lineCap = "round";
    for (let q = 0; q < NB; q++) { const a = (q+0.5)/NB;
      ctx.strokeStyle = `rgba(${col[0]},${col[1]},${col[2]},${(0.3+0.65*a).toFixed(3)})`; ctx.lineWidth = lw; ctx.stroke(cp[q]);
      ctx.strokeStyle = `rgba(${foamCol[0]},${foamCol[1]},${foamCol[2]},${(0.35+0.65*a).toFixed(3)})`; ctx.lineWidth = lw*2.4; ctx.stroke(fp[q]); }
    return true;
  }

  window.FalSwell = {
    init(coast, bathy){ G = build(coast, bathy); resolve(); return G; },
    setSwell(sw){ const changed = !SW || !sw || Math.round(SW.dir/5) !== Math.round(sw.dir/5) || Math.round(SW.tp*2) !== Math.round(sw.tp*2);
      SW = sw; if (changed) resolve(); else shade(); },
    ready(){ return !!(G && F && F.crest); },
    grid(){ return G; },
    draw, CELL
  };
})();
