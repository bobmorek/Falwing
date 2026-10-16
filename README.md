# FAL Wing

Wing foil conditions around Falmouth: Pendennis Point, Gylly, Swanpool, Mylor, Helford mouth and Grebe.

Single-page app (`index.html`), no build step. Live data from Open-Meteo (wind, swell, modelled tide) and the RNAS Culdrose METAR via aviationweather.gov. Measured waves from the Porthleven and Looe Bay buoys (Regional Coastal Monitoring Programme, OGL), which need a free API key in `CCO_KEY`. Map tiles from CARTO/OpenStreetMap with OpenSeaMap seamarks.

Spot settings (shore facing, swell exposure, tidal flag, notes) live in the `SPOTS` list at the top of the script.

Not a safety service. Tide is modelled, not an official prediction.
