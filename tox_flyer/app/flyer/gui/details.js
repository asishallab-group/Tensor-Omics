"use strict";

import {
  createSingleDetailsTable,
  createMasterTable,
  createRowSelector,
  createElement,
  createTableUI,
  createButton,
  show
} from "./dom.js";
import {
  createInputForHeaderData,
  createCustomizationTable,
  applyChanges
} from "./customization.js";

export function setupDetailView() {
  createPickedDetailsDialog();
}

function createTableUIWithCustomizationButton(tableComponents) {
  const button = createButton({
    innerText: "Customize",
    disabled: true,
  });
  button.addEventListener("click", function () {
    show(createTableUI(createCustomizationTable(tableComponents.table), { beforePageSwitch: applyChanges }), applyChanges);
  })

  return createTableUI(tableComponents, {
    bottomChildren: [button],
    beforePageSwitch: applyChanges,
    afterRowSelection(evt) {
      button.disabled = !tableComponents.table.TOX_allSelectorState && tableComponents.table.TOX_selectedRows.size === 0;
    }
  });
}

function multiViewMenuListener({ target }) {
  for (const content of target.closest("section").querySelector("section").children) {
    if (!content.hidden) {
      applyChanges(content);
    }
    content.hidden = content.getAttribute("menu-name") !== target.textContent;
  }
  for (const tab of target.closest("section").querySelector("menu").children) {
    tab.classList.toggle("active", tab.textContent === target.textContent);
  }
}

function createMultiView(contents) {
  const menu = createElement("menu", {classes: ["tabs"]});
  const contentBody = createElement("section");
  for (const [menuName, content] of Object.entries(contents)) {
    const menuItem = createElement("li", {
      textContent: menuName,
      classes: ["clickable"]
    });
    menuItem.addEventListener("click", multiViewMenuListener)
    menu.appendChild(menuItem);

    content.hidden = true;
    content.setAttribute("menu-name", menuName);
    contentBody.appendChild(content);
  }

  if (contentBody.firstChild !== undefined) {
    contentBody.firstChild.hidden = false;
    menu.firstChild.classList.add("active");
  }
  return createElement("section", {
    children: [
      menu,
      contentBody
    ]
  });
}

function createPickedDetailsDialog() {
  if (!document.getElementById("pickedDetails")) {
    {
      const tableHeaders = getDetailsTableDataMap("Gene", "ShiftVector", "Centroid");
      const menuNames = {Gene: "Genes", ShiftVector: "Shift Vectors", Centroid: "Centroids"};
      const multiViewArg = {};

      for (const [pickType, headers] of Object.entries(tableHeaders)) {
        const { table } = createMasterTable({ elements: [] }, null, headers);
        table.id = pickType + "DetailsTable";

        let tableUI;
        if (pickType === "Gene") {
          tableUI = createTableUIWithCustomizationButton({ table });
        } else {
          tableUI = createTableUI({ table }, { beforePageSwitch: applyChanges });
        }

        multiViewArg[menuNames[pickType]] = tableUI;
      }
      const pickedDetails = createElement("dialog", {
        id: "pickedDetails",
        children: [createMultiView(multiViewArg)]
      });
      pickedDetails.addEventListener("close", () => { applyChanges(pickedDetails); config.update() });
      document.getElementById("UI")?.appendChild(pickedDetails);
    }
  }
}

export function appendDetailRow({ family, gene, type }) {
  const id = `${family}.${gene}.${type}`;
  if (document.getElementById(id) === null) {
    let data;
    if (type === "Centroid") {
      const familyData = dataHandler.getFamilyData(family, ...dataHandler.tissues);
      data = { coordinates: familyData.centroid, family: familyData.family };
    } else {
      data = dataHandler.getGeneData(family, gene);
    }
    const table = document.getElementById(type + "DetailsTable");
    table.TOX_elements.push({ family, gene });
    const row = createElement("tr", { id, "tox-family": family, "tox-gene": gene });
    row.appendChild(createRowSelector(false, family, gene, type));

    const dataMap = getDetailsTableDataMap(type);
    for (const cell of dataMap[type]) {
      const td = createElement("td", {
        children: [cell.data(data, family, gene)]
      });
      row.appendChild(td);
    }
    table.tBodies[0].appendChild(row);
    table.querySelector("#selectAll").checked = false;

    return true;
  }

  return false;
}

export function removeDetailRow({ family, gene, type }) {
  const row = document.getElementById(`${family}.${gene}.${type}`);
  if (row !== null) {
    row.querySelector(".row-selector").click();
    const table = row.closest("table");
    const elementIdx = table.TOX_elements.findIndex(element => (family === element.family) && (gene === element.gene));
    table.TOX_elements.splice(elementIdx, 1);
    row.remove();
    return true;
  }

  return false;
}

function createCellLinkElement(value, linkContent) {
  const a = createElement("a", {
    textContent: value,
    classes: ["clickable"]
  });
  a.addEventListener("click", evt => {
    evt.preventDefault();
    applyChanges(evt.target.closest("table"));
    show(linkContent());
  });
  return a;
}

function getDetailsTableDataMap(...types) {
  const headers = {};

  const tissues = dataHandler.tissues;
  const links = {
    family: {
      title: "Family",
      data(geneData, familyIdx, geneIdx) {
        return createCellLinkElement(geneData.family, () => createTableUI(createSingleDetailsTable(geneData, familyIdx, geneIdx, headers.Family), {}));
      }
    },
    gene: {
      title: "Gene",
      data(geneData, familyIdx, geneIdx) {
        return createCellLinkElement(geneData.genes, () => createTableUI(createSingleDetailsTable(geneData, familyIdx, geneIdx, [
          { title: "Identifier", data() {return geneData.genes } },
          links.family,
          { title: "Type", data(geneData) {return geneData.is_outlier ? "Outlier" : "Inlier"} },
          { title: "Species", data() {return geneData.species} },
          { title: "Description", data() {return "..."} },
          ...tissues.map((tissue, i) => {
            return { title: tissue, data() {return geneData.coordinates[i]} }
          })
        ]), {}));
      }
    }
  }
  function tissueRelatedHeader(key) {
    return {
      get title() {
        return createElement("span", {
          classes: [key],
          innerText: config.get(key)
        });
      },
      data(geneData) {
        return createElement("span", {
          classes: [key, "data"],
          innerText: geneData.coordinates[tissues.indexOf(config.get(key))].toFixed(3)
        });
      }
    };
  }
  const tissueRelated = [
    tissueRelatedHeader("tissueX"),
    tissueRelatedHeader("tissueY"),
    tissueRelatedHeader("tissueZ")
  ];

  headers.Gene = [
    links.gene,
    links.family,
    { title: "Type", data(geneData) {return geneData.is_outlier ? "Outlier" : "Inlier"} },
    ...tissueRelated
  ];

  headers.Family = [
    { title: "Identifier", data(geneData) {return geneData.family} },
    {
      title: "Genes",
      data(geneData, familyIdx) {
        return createCellLinkElement(`Inspect ${dataHandler.getGeneCount(familyIdx)} members`, () => createTableUIWithCustomizationButton(createMasterTable(
          { elements: dataHandler.genes(familyIdx).map(gene => ({family: familyIdx, gene})) },
          (family, gene) => dataHandler.getGeneData(family, gene),
          headers.Gene
        )))
      }
    },
    { title: "Description", data() {return "..."}},
  ];

  const dataMap = {};

  for (const type of types) {
    switch (type) {
      case "Gene": {
        dataMap.Gene = headers.Gene;
        break;
      }
      case "ShiftVector": {
        dataMap.ShiftVector = [
          links.gene,
          links.family,
          { title: "Visibility", data: createInputForHeaderData("ShiftVector", "boolean") }
        ];
        break;
      }
      case "Centroid": {
        dataMap.Centroid = [
          links.family,
          { title: "Visibility", data: createInputForHeaderData("Centroid", "boolean") },
          ...tissueRelated
        ];
        break;
      }
      case "Family": {
        dataMap.Family = headers.Family;
      }
    }
  }
  return dataMap;
}
