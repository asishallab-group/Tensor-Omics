d3.json("json/new_multi.json").then(function (data) {

// const svg = d3.select("#figure_tissue_change");
    const margin = { top: 40, right: 20, bottom: 100, left: 50 },
      width = 1200 - margin.left - margin.right,
      height = 600 - margin.top - margin.bottom;

const svg = d3.select("#figure_tissue_change")
  .append("g")
  .attr("transform", `translate(${margin.left},${margin.top})`);

  const tissues = data.tissues;
  const orthogroups = Object.entries(data.orthogroups);
  const selected_genes = data.selected_genes || [];
  const just_outliers = data.just_outliers || false;

  let filtered = [];

  orthogroups.forEach(([og_name, og]) => {
    const vectors = og.unit_vectors;
    const types = og.gene_type || [];
    const outliers = og.outliers || [];
    const gene_names = og.genes || [];

    vectors.forEach((uv, i) => {
      const is_out = outliers[i];
      const gene_id = gene_names[i] || `${og_name}_${i}`;

      const include = (selected_genes.length > 0 && selected_genes.includes(gene_id)) ||
                      (selected_genes.length === 0 && (!just_outliers || is_out));

      if (include) {
        const type = (gene_names.length && types.length && gene_names[i] && types[i]) ? types[i] : "unknown";
        uv.forEach((val, j) => {
          filtered.push({
            tissue: tissues[j],
            value: val,
            gene_type: type
          });
        });
      }
    });
  });

  const grouped = d3.groups(filtered, d => d.tissue, d => d.gene_type);
  const tissueTypes = Array.from(new Set(filtered.map(d => d.tissue)));
  const geneTypes = Array.from(new Set(filtered.map(d => d.gene_type)));

  const x0 = d3.scaleBand()
    .domain(tissueTypes)
    .range([0, width])
    .paddingInner(0.1);

  const x1 = d3.scaleBand()
    .domain(geneTypes)
    .range([0, x0.bandwidth()])
    .padding(0.2);

  const allValues = filtered.map(d => d.value);
  const y = d3.scaleLinear()
    .domain([d3.min(allValues), d3.max(allValues)])
    .nice()
    .range([height, 0]);

  svg.append("g")
    .call(d3.axisLeft(y));

  svg.append("g")
    .attr("transform", `translate(0,${height})`)
    .call(d3.axisBottom(x0))
    .selectAll("text")
    .attr("transform", "rotate(45)")
    .style("text-anchor", "start");

  const boxplot = svg.selectAll("g.boxgroup")
    .data(grouped)
    .enter()
    .append("g")
    .attr("transform", d => `translate(${x0(d[0])},0)`);

  const colorScale = d3.scaleOrdinal()
    .domain(["ortholog", "paralog", "unknown"])
    .range(["green", "orange", "gray"]);

  boxplot.selectAll("g.subbox")
    .data(d => d[1])
    .enter()
    .append("g")
    .attr("transform", d => `translate(${x1(d[0])},0)`)
    .each(function(d) {
      const values = d[1].map(v => v.value).sort(d3.ascending);
      const q1 = d3.quantile(values, 0.25);
      const q2 = d3.quantile(values, 0.5);
      const q3 = d3.quantile(values, 0.75);
      const iqr = q3 - q1;
      const min = Math.max(d3.min(values), q1 - 1.5 * iqr);
      const max = Math.min(d3.max(values), q3 + 1.5 * iqr);
      const g = d3.select(this);

      // Draw box
      g.append("rect")
        .attr("y", y(q3))
        .attr("height", y(q1) - y(q3))
        .attr("width", x1.bandwidth())
        .attr("fill", colorScale(d[0]))
        .attr("stroke", "black");

      // Draw median line
      g.append("line")
        .attr("x1", 0)
        .attr("x2", x1.bandwidth())
        .attr("y1", y(q2))
        .attr("y2", y(q2))
        .attr("stroke", "black")
        .attr("stroke-width", 1.5);

      // Draw min/max whiskers
      g.append("line")
        .attr("x1", x1.bandwidth() / 2)
        .attr("x2", x1.bandwidth() / 2)
        .attr("y1", y(min))
        .attr("y2", y(max))
        .attr("stroke", "black");

      // Horizontal lines on whiskers
      g.append("line")
        .attr("x1", 0)
        .attr("x2", x1.bandwidth())
        .attr("y1", y(min))
        .attr("y2", y(min))
        .attr("stroke", "black");
      g.append("line")
        .attr("x1", 0)
        .attr("x2", x1.bandwidth())
        .attr("y1", y(max))
        .attr("y2", y(max))
        .attr("stroke", "black");
    });
});

  