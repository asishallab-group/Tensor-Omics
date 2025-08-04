/**
 * JSON Config Schema for generixXYPlot
 *
 * Required fields depend on plotType:
 *
 * plotType: 'xy' | 'histogram' | 'bar' (default: 'xy')
 *
 * For plotType: 'xy' (points/lines)
 *   - xValues: [number|string]   // required
 *   - yValues: [number]         // required
 *   - xField: string            // optional (default: 'x')
 *   - yField: string            // optional (default: 'y')
 *   - showPoints: boolean       // optional (default: true)
 *   - connectPoints: boolean    // optional (default: false)
 *   - lineStyle: 'straight' | 'smooth' // optional (default: 'straight')
 *   - lineColor: string         // optional (default: '#0074D9')
 *   - pointColor: string        // optional (default: '#FF4136')
 *   - highlightPoints: {x:[], y:[], color, label} // optional
 *   - highlightIndices: [int]   // optional
 *   - highlightXValues: [number]// optional
 *   - highlightColor: string    // optional (default: '#FFD700')
 *   - xLabel: string            // optional
 *   - yLabel: string            // optional
 *   - yMin, yMax: number        // optional (force y axis limits)
 *   - width, height: number     // optional (default: 600x400)
 *   - margin: {top,right,bottom,left} // optional
 *
 * For plotType: 'histogram'
 *   - yValues: [number]         // required
 *   - yField: string            // optional (default: 'y')
 *   - bins: number              // optional (default: 20)
 *       (Note: The actual number of bars may differ from 'bins' due to D3's binning algorithm and data distribution.)
 *   - barColor: string          // optional (default: '#0074D9')
 *   - xLabel, yLabel: string    // optional
 *   - width, height, margin     // optional
 *
 * For plotType: 'bar'
 *   - xValues: [string|number]  // required
 *   - yValues: [number]         // required
 *   - xField: string            // optional (default: 'x')
 *   - yField: string            // optional (default: 'y')
 *   - barColor: string          // optional (default: '#0074D9')
 *   - xLabel, yLabel: string    // optional
 *   - width, height, margin     // optional
 *
 * Example:
 * {
 *   "plotType": "xy",
 *   "xValues": [0,1,2],
 *   "yValues": [10,5,30],
 *   "showPoints": true,
 *   "connectPoints": true,
 *   "lineStyle": "smooth",
 *   "highlightPoints": {"x": [1], "y": [20], "color": "#FFD700", "label": "Ref"},
 *   "xLabel": "Time (s)",
 *   "yLabel": "Value"
 * }
 *
 * {
 *   "plotType": "histogram",
 *   "yValues": [1,2,2,3,3,3,4],
 *   "bins": 3,
 *   "xLabel": "Value",
 *   "yLabel": "Count"
 * }
 *
 * {
 *   "plotType": "bar",
 *   "xValues": ["A","B","C"],
 *   "yValues": [5,8,2],
 *   "xLabel": "Category",
 *   "yLabel": "Value"
 * }
 */

// Load the JSON config and data, and draw the plot automatically
document.addEventListener('DOMContentLoaded', function() {
  d3.json('json/test.json').then(config => {
    generixXYPlot('#xy_plot', config);
  });
});

/**
 * Generic D3.js XY Plot Component
 * Usage:
 *   generixXYPlot(container, data, config)
 *
 * @param {HTMLElement|string} container - DOM element or selector
 * @param {Array<Object>} data - Array of objects with x/y fields
 * @param {Object} config - JSON config (see issue description)
 */

function generixXYPlot(container, dataOrConfig, config) {
  // Allows call: generixXYPlot(container, configWithValues) or generixXYPlot(container, data, config)
  let data, cfg;
  const defaultColors = {
    lineColor: '#0074D9',
    pointColor: '#FF4136',
    highlightColor: '#FFD700'
  };
  if (Array.isArray(dataOrConfig)) {
    // Classic mode: data, config
    data = dataOrConfig;
    cfg = Object.assign({
      plotType: 'xy', // 'xy' (default), 'histogram', 'bar'
      xField: 'x',
      yField: 'y',
      showPoints: true,
      connectPoints: false,
      lineStyle: 'straight',
      width: 600,
      height: 400,
      margin: {top: 40, right: 30, bottom: 50, left: 60},
      xLabel: null,
      yLabel: null,
      ...defaultColors
    }, config);
  } else {
    // Only config with xValues/yValues or yValues only (histogram)
    cfg = Object.assign({
      plotType: 'xy', // 'xy' (default), 'histogram', 'bar'
      xField: 'x',
      yField: 'y',
      showPoints: true,
      connectPoints: false,
      lineStyle: 'straight',
      width: 600,
      height: 400,
      margin: {top: 40, right: 30, bottom: 50, left: 60},
      xLabel: null,
      yLabel: null,
      ...defaultColors
    }, dataOrConfig);
    if (cfg.plotType === 'histogram') {
      // Only yValues needed for histogram
      if (Array.isArray(cfg.yValues)) {
        data = cfg.yValues.map(y => ({ [cfg.yField]: y }));
      } else {
        throw new Error('Config for histogram must include yValues as array');
      }
    } else if (cfg.plotType === 'bar') {
      // xValues and yValues needed for bar
      if (Array.isArray(cfg.xValues) && Array.isArray(cfg.yValues)) {
        data = cfg.xValues.map((x, i) => ({
          [cfg.xField]: x,
          [cfg.yField]: cfg.yValues[i]
        }));
      } else {
        throw new Error('Config for bar must include xValues and yValues as arrays');
      }
    } else {
      // Default: xy plot needs both xValues and yValues
      if (Array.isArray(cfg.xValues) && Array.isArray(cfg.yValues)) {
        data = cfg.xValues.map((x, i) => ({
          [cfg.xField]: x,
          [cfg.yField]: cfg.yValues[i]
        }));
      } else {
        throw new Error('Config for xy plot must include xValues and yValues as arrays');
      }
    }
  }

  // Responsive: clear previous SVG
  const sel = (typeof container === 'string') ? d3.select(container) : d3.select(container);
  sel.select('svg').remove();

  const width = cfg.width - cfg.margin.left - cfg.margin.right;
  const height = cfg.height - cfg.margin.top - cfg.margin.bottom;

  // SVG setup
  const svg = sel.append('svg')
    .attr('width', cfg.width)
    .attr('height', cfg.height)
    .attr('viewBox', `0 0 ${cfg.width} ${cfg.height}`)
    .attr('preserveAspectRatio', 'xMidYMid meet');

  const g = svg.append('g')
    .attr('transform', `translate(${cfg.margin.left},${cfg.margin.top})`);

  // Plot type switch: 'xy' (default), 'histogram', 'bar'
  const plotType = cfg.plotType || 'xy';

  if (plotType === 'histogram') {
    // Histogram: use yField as the variable to bin
    const values = data.map(d => +d[cfg.yField]);
    const bins = d3.bin().thresholds(cfg.bins || 20)(values);
    const x = d3.scaleLinear()
      .domain([bins[0].x0, bins[bins.length - 1].x1])
      .range([0, width]);
    const y = d3.scaleLinear()
      .domain([0, d3.max(bins, d => d.length)])
      .range([height, 0]);

    // Axes
    g.append('g')
      .attr('class', 'x-axis')
      .attr('transform', `translate(0,${height})`)
      .call(d3.axisBottom(x));
    g.append('g')
      .attr('class', 'y-axis')
      .call(d3.axisLeft(y));

    // Bars
    g.selectAll('.bar')
      .data(bins)
      .enter().append('rect')
      .attr('class', 'bar')
      .attr('x', d => x(d.x0))
      .attr('y', d => y(d.length))
      .attr('width', d => Math.max(0, x(d.x1) - x(d.x0) - 1))
      .attr('height', d => height - y(d.length))
      .attr('fill', cfg.barColor || '#0074D9');

    // Axis labels
    g.append('text')
      .attr('class', 'x-label')
      .attr('x', width/2)
      .attr('y', height + 40)
      .attr('text-anchor', 'middle')
      .text(cfg.xLabel || cfg.yField + ' value');
    g.append('text')
      .attr('class', 'y-label')
      .attr('x', -height/2)
      .attr('y', -45)
      .attr('transform', 'rotate(-90)')
      .attr('text-anchor', 'middle')
      .text('Count');
    return;
  }

  if (plotType === 'bar') {
    // Bar chart: xField as categories, yField as value
    const x = d3.scaleBand()
      .domain(data.map(d => d[cfg.xField]))
      .range([0, width])
      .padding(0.1);
    const y = d3.scaleLinear()
      .domain([0, d3.max(data, d => +d[cfg.yField])])
      .range([height, 0]);

    // Axes
    g.append('g')
      .attr('class', 'x-axis')
      .attr('transform', `translate(0,${height})`)
      .call(d3.axisBottom(x));
    g.append('g')
      .attr('class', 'y-axis')
      .call(d3.axisLeft(y));

    // Bars
    g.selectAll('.bar')
      .data(data)
      .enter().append('rect')
      .attr('class', 'bar')
      .attr('x', d => x(d[cfg.xField]))
      .attr('y', d => y(+d[cfg.yField]))
      .attr('width', x.bandwidth())
      .attr('height', d => height - y(+d[cfg.yField]))
      .attr('fill', cfg.barColor || '#0074D9');

    // Axis labels
    g.append('text')
      .attr('class', 'x-label')
      .attr('x', width/2)
      .attr('y', height + 40)
      .attr('text-anchor', 'middle')
      .text(cfg.xLabel || cfg.xField);
    g.append('text')
      .attr('class', 'y-label')
      .attr('x', -height/2)
      .attr('y', -45)
      .attr('transform', 'rotate(-90)')
      .attr('text-anchor', 'middle')
      .text(cfg.yLabel || cfg.yField);
    return;
  }

  // Default: XY plot (points/lines)
  // Scales
  const x = d3.scaleLinear()
    .domain(d3.extent(data, d => +d[cfg.xField]))
    .range([0, width]);

  // Y domain: min/max of yValues and highlightPoints.y
  let yVals = data.map(d => +d[cfg.yField]);
  if (cfg.highlightPoints && Array.isArray(cfg.highlightPoints.y)) {
    yVals = yVals.concat(cfg.highlightPoints.y.map(Number));
  }
  const yMin = (typeof cfg.yMin === 'number') ? cfg.yMin : d3.min(yVals);
  const yMax = (typeof cfg.yMax === 'number') ? cfg.yMax : d3.max(yVals);
  const y = d3.scaleLinear()
    .domain([yMin, yMax])
    .range([height, 0]);

  // Axes
  g.append('g')
    .attr('class', 'x-axis')
    .attr('transform', `translate(0,${height})`)
    .call(d3.axisBottom(x));
  g.append('g')
    .attr('class', 'y-axis')
    .call(d3.axisLeft(y));

  // Axis labels
  g.append('text')
    .attr('class', 'x-label')
    .attr('x', width/2)
    .attr('y', height + 40)
    .attr('text-anchor', 'middle')
    .text(cfg.xLabel || cfg.xField);
  g.append('text')
    .attr('class', 'y-label')
    .attr('x', -height/2)
    .attr('y', -45)
    .attr('transform', 'rotate(-90)')
    .attr('text-anchor', 'middle')
    .text(cfg.yLabel || cfg.yField);

  // Line generator
  let curveType = d3.curveLinear;
  if (cfg.lineStyle === 'smooth') curveType = d3.curveMonotoneX;
  const line = d3.line()
    .x(d => x(+d[cfg.xField]))
    .y(d => y(+d[cfg.yField]))
    .curve(curveType);

  // Sort data by x for line
  const sortedData = [...data].sort((a, b) => +a[cfg.xField] - +b[cfg.xField]);

  // Draw line if needed
  if (cfg.connectPoints) {
    g.append('path')
      .datum(sortedData)
      .attr('fill', 'none')
      .attr('stroke', cfg.lineColor)
      .attr('stroke-width', 2)
      .attr('d', line);
  }

  // Draw main points, with optional highlighting by index or x value
  if (cfg.showPoints) {
    const highlight = Array.isArray(cfg.highlightIndices) ? cfg.highlightIndices : [];
    const highlightX = Array.isArray(cfg.highlightXValues) ? cfg.highlightXValues : [];
    g.selectAll('.point')
      .data(data)
      .enter().append('circle')
      .attr('class', (d, i) => (highlight.includes(i) || highlightX.includes(+d[cfg.xField])) ? 'point highlight' : 'point')
      .attr('cx', d => x(+d[cfg.xField]))
      .attr('cy', d => y(+d[cfg.yField]))
      .attr('r', 4)
      .attr('fill', (d, i) => (highlight.includes(i) || highlightX.includes(+d[cfg.xField])) ? cfg.highlightColor : cfg.pointColor)
      .attr('stroke', '#fff')
      .attr('stroke-width', (d, i) => (highlight.includes(i) || highlightX.includes(+d[cfg.xField])) ? 2.5 : 1);
  }

  // Draw highlightPoints if present
  if (cfg.highlightPoints && Array.isArray(cfg.highlightPoints.x) && Array.isArray(cfg.highlightPoints.y)) {
    const hp = cfg.highlightPoints;
    g.selectAll('.highlight-point')
      .data(hp.x.map((xv, i) => ({ x: xv, y: hp.y[i] })))
      .enter().append('circle')
      .attr('class', 'highlight-point')
      .attr('cx', d => x(+d.x))
      .attr('cy', d => y(+d.y))
      .attr('r', 5)
      .attr('fill', hp.color || '#FFD700')
      .attr('stroke', '#222')
      .attr('stroke-width', 2);
  }
}

// Example usage:
// generixXYPlot('#plot', data, config);
// To update, just call again with new data/config (it will clear previous SVG)
// Points can be highlighted by passing highlightIndices: [i1, i2, ...] in config
