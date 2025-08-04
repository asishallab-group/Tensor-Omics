d3.json("json/figures_mult.json").then(function (data) {
  const svg = d3.select("#figure_b_svg");
  const width = 1000;
  const height = 500;
  const margin = 50;
  const centerX = width / 2;
  const centerY = height / 2;

  // Color scale for each orthogroup/family
  const colorScale = d3.scaleOrdinal(d3.schemeCategory10);
  const just_outliers = data.just_outliers || false; 

  // Extract tissues and orthogroups
  const tissues = data.tissues; // Extract the tissue list
  const orthogroups = Object.keys(data.orthogroups || {});

  // Scale to convert angle to radius for visualization
  const radiusScale = d3.scaleLinear()
    .domain([0, Math.PI])
    .range([50, Math.min(width, height) / 2 - margin]);

  // Define diagonal vector and its magnitude for angle calculation
  const dimensions = tissues.length;
  const diagonalVector = Array(dimensions).fill(1 / Math.sqrt(dimensions));
  const magnitudeDiagonal = Math.sqrt(diagonalVector.reduce((sum, val) => sum + val ** 2, 0));

  // Function to compute angle between a vector and the diagonal
  const getAngle = (vec) => {
    const dot = vec.reduce((sum, val, i) => sum + val * diagonalVector[i], 0);
    const magVec = Math.sqrt(vec.reduce((sum, val) => sum + val ** 2, 0));
    return Math.acos(Math.min(1, Math.max(-1, dot / (magVec * magnitudeDiagonal))));
  };

  // Convert polar coordinates to cartesian
  const polarToCartesian = (r, angle) => ({
    x: centerX + r * Math.sin(angle),
    y: centerY - r * Math.cos(angle)
  });

  let maxOverflow = 0;

  // Define a mapping of specific families to their colors
  const familyColors = {
    "OG0000017": "red",
    "OG0000034": "green",
    "OG0000306": "gray",
    "OG0000103": "cyan",
    "OG0000141": "magenta"
  };
  // Precompute all needed info per orthogroup
  const allGroups = orthogroups.map((orthogroup, familyIndex) => {

    const { vectors, centroid, outliers } = data.orthogroups[orthogroup];
    const angles = vectors.map(getAngle);
    const centroidAngle = getAngle(centroid);

    const positions = angles.map((angle, index) => {
      const radius = radiusScale(angle) + index * 15;
      const { x, y } = polarToCartesian(radius, angle);
      if (y < margin) maxOverflow = Math.max(maxOverflow, margin - y);
      return { angle, radius, x, y, isOutlier: outliers[index], tissue: tissues[index] }; // Include tissue info
    });

    const centroidRadius = radiusScale(centroidAngle) + angles.length * 15;
    const { x: cx, y: cy } = polarToCartesian(centroidRadius, centroidAngle);
    if (cy < margin) maxOverflow = Math.max(maxOverflow, margin - cy+200);

    // Use the specific color for the family if defined, otherwise use the default color scale
    const familyColor = familyColors[orthogroup] || colorScale(familyIndex);

    return {
      orthogroup,
      familyIndex,
      familyColor,
      // familyColor: colorScale(familyIndex),
      positions,
      centroidAngle,
      centroidRadius,
      centroidX: cx,
      centroidY: cy
    };
  });

  // Main plotting group with vertical shift to avoid top overflow
  const plotGroup = svg.append("g").attr("transform", `translate(0, ${maxOverflow})`);
  const totalHeight = height + maxOverflow;
d3.select("#figure_b_svg")
  .attr("height", totalHeight);

  const defs = svg.append("defs"); // Marker definitions for arrows

  // Draw vectors, arcs, centroid and markers
  allGroups.forEach(group => {
    const g = plotGroup.append("g").attr("class", `family-${group.orthogroup}`);
    const { familyColor, positions, centroidAngle, centroidRadius, centroidX, centroidY, familyIndex } = group;

    // Arrow marker for centroid
    defs.append("marker")
      .attr("id", `triangle-${familyIndex}`)
      .attr("viewBox", "0 0 10 10")
      .attr("refX", 5).attr("refY", 5)
      .attr("markerWidth", 6).attr("markerHeight", 6)
      .attr("orient", "auto")
      .append("path")
      .attr("d", "M 0 0 L 10 5 L 0 10 Z")
      .attr("fill", familyColor);

    // Centroid arc
    g.append("path")
      .attr("d", d3.arc()
        .innerRadius(centroidRadius - 10)
        .outerRadius(centroidRadius)
        .startAngle(0)
        .endAngle(centroidAngle)())
      .attr("transform", `translate(${centerX}, ${centerY})`)
      .attr("fill", "#66A61E")
      .attr("opacity", 0.5);

    // Centroid line
    g.append("line")
      .attr("x1", centerX)
      .attr("y1", centerY)
      .attr("x2", centroidX)
      .attr("y2", centroidY)
      .attr("stroke", familyColor)
      .attr("stroke-width", 5)
      .attr("marker-end", `url(#triangle-${familyIndex})`);

    // Vectors and arcs
    positions.forEach(({ angle, radius, x, y, isOutlier, tissue }) => {
      if (!just_outliers) {
        g.append("line")
          .attr("x1", centerX)
          .attr("y1", centerY)
          .attr("x2", x)
          .attr("y2", y)
          .attr("stroke", familyColor)
          .attr("stroke-width", 2);

        g.append("path")
          .attr("d", d3.arc()
            .innerRadius(radius - 3)
            .outerRadius(radius)
            .startAngle(0)
            .endAngle(angle)())
          .attr("transform", `translate(${centerX}, ${centerY})`)
          .attr("fill", familyColor)
          .attr("opacity", 0.5);
      }
      if (isOutlier) {
        g.append("line")
          .attr("x1", centerX)
          .attr("y1", centerY)
          .attr("x2", x)
          .attr("y2", y)
          .attr("stroke", familyColor)
          .attr("stroke-width", 2);
          
        g.append("path")
        .attr("d", d3.arc()
          .innerRadius(radius - 3)
          .outerRadius(radius)
          .startAngle(0)
          .endAngle(angle)())
        .attr("transform", `translate(${centerX}, ${centerY})`)
        .attr("fill", familyColor)
        .attr("opacity", 0.5);

        g.append("circle")
          .attr("cx", x)
          .attr("cy", y)
          .attr("r", 5)
          .attr("fill", familyColor);

      }
    });
  });

  // Angle labels group
  const labelsGroup = svg.append("g")
    .attr("class", "labels-group")
    .attr("transform", `translate(0, ${maxOverflow})`);

  // Add text labels for angles and centroid
  allGroups.forEach(group => {
    const { positions, centroidAngle, centroidRadius } = group;

    positions.forEach(({ angle, radius, isOutlier }) => {
      const labelRadius = radius - 5;
      const { x, y } = polarToCartesian(labelRadius, angle / 2);
      if(isOutlier && just_outliers || !just_outliers) {
        labelsGroup.append("text")
          .attr("x", x)
          .attr("y", y)
          .text(`${(angle * 180 / Math.PI).toFixed(2)}°`)
          .attr("fill", "black")
          .style("font-size", "12px")
          .style("text-anchor", "middle");
      }
      // labelsGroup.append("text")
      //   .attr("x", x)
      //   .attr("y", y)
      //   .text(`${(angle * 180 / Math.PI).toFixed(2)}°`)
      //   .attr("fill", "black")
      //   .style("font-size", "12px")
      //   .style("text-anchor", "middle");
    });

    const labelRadius = centroidRadius - 10;
    const { x, y } = polarToCartesian(labelRadius, centroidAngle / 2);
    labelsGroup.append("text")
      .attr("x", x)
      .attr("y", y)
      .text(`${(centroidAngle * 180 / Math.PI).toFixed(2)}°`)
      .attr("fill", "black")
      .style("font-size", "12px")
      .style("text-anchor", "middle");
  });

  // Draw diagonal axis (reference vector)
  plotGroup.append("line")
    .attr("x1", centerX)
    .attr("y1", Math.min(margin - maxOverflow, 0))
    .attr("x2", centerX)
    .attr("y2", height - margin)
    .attr("stroke", "deeppink")
    .attr("stroke-width", 5);
});