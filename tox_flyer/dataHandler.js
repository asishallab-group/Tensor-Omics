"use strict";

import data from "./sampleData.json" with { type: "json" };

export const handler = {
  get families() {
    return Object.keys(data);
  },
  get tissues() {
    return Object.keys(
      data[this.families[0]]?.[tissues] ?? {}
    )
  },
  getCentroid(family, ...tissues) {
    return tissues.map((tissue) => {
      return data[family]?.["centroid"][tissue];
    })
  },
  iterGenes: function* (family, ...tissues) {
    const familyData = data[family];
    if (familyData !== undefined) {
      const allGeneData = Object.entries(familyData);
      allGeneData.pop(); // remove centroid data
      const allTissueData = allGeneData.pop()[1];
      for (const i in familyData.genes) {
        const singleGeneData = {}
        for (const [key, values] of allGeneData) {
          singleGeneData[key] = values[i];
        }
        singleGeneData["coordinates"] = tissues.map((tissue) => {
          return allTissueData[tissue]?.[i];
        })
        yield singleGeneData;
      }
    }
  }
}
Object.freeze(handler);
