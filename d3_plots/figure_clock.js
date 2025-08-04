d3.json("json/figures_mult.json").then(function (data) {
    const svg = d3.select("#figure_clock_svg");
    const width = 1000;
    const height = 1000;
    const margin = 50;
    const centerX = width / 2;
    const centerY = height / 2;
    // Variable para controlar si solo se dibujan los outliers
    const just_outliers = data.just_outliers || false; 

    // Extract tissues and orthogroups
    const tissues = data.tissues || [];
    const orthogroups = Object.keys(data.orthogroups || {});
    const numTissues = tissues.length;
  
    // Define a color scale for the families
    const colorScale = d3.scaleOrdinal(d3.schemeCategory10);
  
    // Calculate angles for each tissue
    const angleStep = 360 / numTissues; // Divide 360 degrees by the number of tissues
  
    // Define a marker for arrows
    svg.append("defs").append("marker")
    .attr("id", "arrowhead")
    .attr("viewBox", "0 0 10 10")
    .attr("refX", 5) // Position of the arrow tip
    .attr("refY", 5)
    .attr("markerWidth", 6)
    .attr("markerHeight", 6)
    .attr("orient", "auto")
    .append("path")
    .attr("d", "M 0 0 L 10 5 L 0 10 Z") // Arrow shape
    .attr("fill", "brown"); // Arrow color

    // Function to convert polar coordinates to cartesian
    const polarToCartesian = (r, angle) => {
      const radians = (angle * Math.PI) / 180; // Convert degrees to radians
      return {
        x: centerX + r * Math.sin(radians),
        y: centerY - r * Math.cos(radians),
      };
    };
  
    // Function to project a vector to the diagonal
    const projectToDiagonal = (vec) => {
      const n = vec.length;
      const diagonal = Array(n).fill(1 / Math.sqrt(n));
      const dotProduct = vec.reduce((sum, val, i) => sum + val * diagonal[i], 0);
      return diagonal.map(d => dotProduct * d);
    };
  
    // Function to project a vector to the plane perpendicular to the diagonal
    const projectToPlane = (vec) => {
      const projectedToDiagonal = projectToDiagonal(vec);
      return vec.map((val, i) => val - projectedToDiagonal[i]); // Subtract the diagonal component
    };
  
    // Function to normalize a vector
    const normalizeVector = (vec) => {
        const magnitude = Math.sqrt(vec.reduce((sum, val) => sum + val ** 2, 0));
        return vec.map(val => val / magnitude);
    };
    
    // Draw clock-like axes for each tissue
    tissues.forEach((tissue, index) => {
      const angle = index * angleStep; // Calculate the angle for this tissue
      const { x, y } = polarToCartesian(height / 2 - margin, angle); // End point of the axis
  
      // Draw the axis line
      svg.append("line")
        .attr("x1", centerX)
        .attr("y1", centerY)
        .attr("x2", x)
        .attr("y2", y)
        .attr("stroke", "gray")
        .attr("stroke-width", 2);
  
      // Add a label for the tissue
      const labelOffset = 20; // Offset for the label
      const { x: labelX, y: labelY } = polarToCartesian(height / 2 - margin + labelOffset, angle);
      svg.append("text")
        .attr("x", labelX)
        .attr("y", labelY)
        .text(tissue)
        .attr("fill", "black")
        .style("font-size", "12px")
        .style("text-anchor", "middle");
    });
  
        // Add a pink point at the center to represent the diagonal
        svg.append("circle")
        .attr("cx", centerX)
        .attr("cy", centerY)
        .attr("r", 15) // Radius of the circle
        .attr("fill", "deeppink");
    
        // Add a label "D" at the center
        svg.append("text")
        .attr("x", centerX)
        .attr("y", centerY + 6) // Slight adjustment to center the text vertically
        .text("D")
        .attr("fill", "white")
        .style("font-size", "20px")
        .style("font-weight", "bold")
        .style("text-anchor", "middle");

    // Find the maximum magnitude among all vectors
    let maxMagnitude = 0;
    orthogroups.forEach((orthogroup) => {
      const { vectors } = data.orthogroups[orthogroup];
      vectors.forEach((vec) => {
        const projectedToPlane = projectToPlane(vec);
        const magnitude = Math.sqrt(projectedToPlane.reduce((sum, val) => sum + val ** 2, 0));
        maxMagnitude = Math.max(maxMagnitude, magnitude);
      });
    });
  
    // Define the maximum radius for the plot
    const maxRadius = height / 2 - margin;
    const familyColors = {
      "OG0000017": "red",
      "OG0000034": "green",
      "OG0000306": "gray",
      "OG0000103": "cyan",
      "OG0000141": "magenta"
    };
  
    // Process each orthogroup and project vectors to the plane
    orthogroups.forEach((orthogroup, familyIndex) => {
      const { vectors, centroid, outliers } = data.orthogroups[orthogroup]; // Include centroid and outliers array
      // const familyColor = colorScale(familyIndex); // Assign a unique color to the family
      const familyColor = familyColors[orthogroup] || colorScale(familyIndex);

  
      // Draw the centroid as a triangle
      const projectedCentroid = projectToPlane(centroid);
      const centroidMagnitude = Math.sqrt(projectedCentroid.reduce((sum, val) => sum + val ** 2, 0));
      const normalizedCentroidMagnitude = (centroidMagnitude / maxMagnitude) * maxRadius;
      const centroidAngle = Math.atan2(projectedCentroid[1], projectedCentroid[0]) * (180 / Math.PI);
      const { x: cx, y: cy } = polarToCartesian(normalizedCentroidMagnitude, centroidAngle);
  
      svg.append("path")
        .attr("d", d3.symbol().type(d3.symbolTriangle).size(150)()) // Triangle symbol
        .attr("transform", `translate(${cx}, ${cy})`)
        .attr("fill", familyColor);
  
      vectors.forEach((vec, vecIndex) => {
        // Project the vector to the plane
        const projectedToPlane = projectToPlane(vec);
  
        // Calculate the magnitude of the projection
        const magnitude = Math.sqrt(projectedToPlane.reduce((sum, val) => sum + val ** 2, 0));
  
        // Normalize the magnitude to fit within the plot
        const normalizedMagnitude = (magnitude / maxMagnitude) * maxRadius;
  
        // Calculate the angle of the projection in the plane
        const angle = Math.atan2(projectedToPlane[1], projectedToPlane[0]) * (180 / Math.PI);
  
        // Convert the normalized magnitude and angle to cartesian coordinates
        const { x, y } = polarToCartesian(normalizedMagnitude, angle);
  
        // Check if the vector is an outlier
        if (outliers[vecIndex]) {
          console.log("outlier");
          // Draw a large circle for outliers
          svg.append("circle")
            .attr("cx", x)
            .attr("cy", y)
            .attr("r", 8) // Larger radius for outliers
            .attr("fill", familyColor);
  
          // Draw a line with an arrow from the centroid to the outlier
          svg.append("line")
            .attr("x1", cx)
            .attr("y1", cy)
            .attr("x2", x)
            .attr("y2", y)
            .attr("stroke", "brown")
            .attr("stroke-width", 2)
            .attr("marker-end", "url(#arrowhead)");

            // Calculate the unit vector
            const unitVector = normalizeVector([x - cx, y - cy]);


            // Scale the unit vector to match the original vector's length
            const scaleFactor = maxRadius / maxMagnitude;
            const unitX = cx + unitVector[0] * 50;
            const unitY = cy + unitVector[1] * 50;

            // Draw the unit vector on top of the original vector
            svg.append("line")
            .attr("x1", cx)
            .attr("y1", cy)
            .attr("x2", unitX)
            .attr("y2", unitY)
            .attr("stroke", "#D4AF37") // Different color for the unit vector
            .attr("stroke-width", 6)
            .attr("stroke-dasharray", "4,2"); // Dashed line for distinction
            
        } else if (!just_outliers) {
            // Draw a small circle for regular points if just_outliers is false
            svg.append("circle")
              .attr("cx", x)
              .attr("cy", y)
              .attr("r", 3) // Smaller radius for regular points
              .attr("fill", familyColor); // Use the family color
          }
      });
    });
    // Add a legend
  const legend = svg.append("g").attr("transform", `translate(${margin}, ${height - margin})`);
  var aux_val=20;
  // Legend for centroid (triangle)
  legend.append("path")
    .attr("d", d3.symbol().type(d3.symbolTriangle).size(150)())
    .attr("transform", `translate(0, ${-height + margin + aux_val})`) // Reemplazado 5 con -height+margin+20
    .attr("fill", "black");
  legend.append("text")
    .attr("x", 20)
    .attr("y", -height+margin+aux_val)
    .text("Centroid")
    .style("font-size", "12px")
    .attr("alignment-baseline", "middle");

  // Legend for outliers (large circle)
  aux_val=aux_val+20;
  legend.append("circle")
    .attr("cx", 0)
    .attr("cy", -height+margin+aux_val)
    .attr("r", 8)
    .attr("fill", "black");
  legend.append("text")
    .attr("x", 20)
    .attr("y", -height+margin+aux_val)
    .text("Outlier")
    .style("font-size", "12px")
    .attr("alignment-baseline", "middle");

    aux_val=aux_val+20;
    if (!just_outliers) {
        // Legend for regular points (small circle)
        legend.append("circle")
            .attr("cx", 0)
            .attr("cy", -height+margin+aux_val)
            .attr("r", 3)
            .attr("fill", "black");
        legend.append("text")
            .attr("x", 20)
            .attr("y", -height+margin+aux_val)
            .text("Regular Gene")
            .style("font-size", "12px")
            .attr("alignment-baseline", "middle");
        aux_val=aux_val+20;
    }

  // Legend for unit vector (dashed line)
  legend.append("line")
    .attr("x1", -10)
    .attr("y1", -height+margin+aux_val)
    .attr("x2", 10)
    .attr("y2", -height+margin+aux_val)
    .attr("stroke", "#D4AF37")
    .attr("stroke-width", 6)
    .attr("stroke-dasharray", "4,2");
  legend.append("text")
    .attr("x", 20)
    .attr("y", -height+margin+aux_val)
    .text("Unit Vector")
    .style("font-size", "12px")
    .attr("alignment-baseline", "middle");
    aux_val=aux_val+20;
    // Legend for shift vector 
    legend.append("line")
    .attr("x1", -10)
    .attr("y1", -height+margin+aux_val)
    .attr("x2", 10)
    .attr("y2", -height+margin+aux_val)
    .attr("stroke", "brown")
    .attr("stroke-width", 2)
    .attr("marker-end", "url(#arrowhead)");
    legend.append("text")
    .attr("x", 20)
    .attr("y", -height+margin+aux_val)
    .text("Shift Vector")
    .style("font-size", "12px")
    .attr("alignment-baseline", "middle");
  });