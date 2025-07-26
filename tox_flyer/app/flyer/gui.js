"use strict";

export function setupGUI() {
  createUIdiv();
  createPickedDetailsDialog();
  switchToDetails("Gene");

  let pickedCount = 0;
  document.addEventListener("pick", evt => {
    if (pickedCount === 0) {
      createDetailButton();
    }
    appendDetailRow(evt.detail);
    pickedCount++;
  });

  document.addEventListener("unpick", evt => {
    pickedCount--;
    removeDetailRow(evt.detail);
    if (pickedCount === 0) {
      removeDetailButton();
    }
  });
}

function appendDetailRow({ family, gene, type }) {
  const id = `${family}.${gene}.${type}`;
  let data;
  if (type === "Centroid") {
    const familyData = dataHandler.getFamilyData(family, ...dataHandler.tissues);
    data = { coordinates: familyData.centroid, family: familyData.family };
  } else {
    data = dataHandler.getGeneData(family, gene);
  }
  const table = document.getElementById(type + "DetailsTable");
  const row = createElement("tr", { id });
  row.appendChild(createRowSelector(false, family, gene, type));

  const dataMap = getDetailsTableDataMap();
  for (const cell of dataMap[type]) {
    const td = createElement("td", {
      children: [cell.data(data, family, gene)]
    });
    row.appendChild(td);
  }
  table.tBodies[0].appendChild(row);
}

function switchToDetails(type) {
  for (const table of document.getElementsByClassName("detailsTable")) {
    table.hidden = table.id !== type + "DetailsTable";
  }
}

function removeDetailRow({ family, gene, type }) {
  document.getElementById(`${family}.${gene}.${type}`)?.remove();
}

function createPickedDetailsDialog() {
  if (!document.getElementById("pickedDetails")) {
    {
      const menu = createElement("menu", { id: "pickedDetailsMenu" });
      const tables = createElement("section", { id: "pickedDetailsTables" });

      const tableHeaders = getDetailsTableDataMap();

      for (const [pickType, headers] of Object.entries(tableHeaders)) {
        const menuItem = createElement("li", {
          id: pickType + "Details",
          textContent: pickType + "s",
          classes: ["clickable"]
        });
        menuItem.addEventListener("click", () => {
          switchToDetails(pickType);
        })
        menu.appendChild(menuItem);

        const table = createMasterTable([], null, headers);
        table.id = pickType + "DetailsTable";
        table.hidden = true;
        table.classList.add("detailsTable");
        tables.appendChild(table);
      }
      const pickedDetails = createElement("dialog", {
        id: "pickedDetails",
        children: [
          menu,
          tables
        ]
      });
      document.getElementById("UI")?.appendChild(pickedDetails);
    }
  }
}

function createRowSelector(selectAll, family, gene, type) {
  const checkbox = createElement("input", { type: "checkbox", classes: ["clickable"] });
  if (selectAll) {
    checkbox.id = "selectAll";
    checkbox.addEventListener("change", evt => {
      const table = evt.target.closest("table");
      const changeEvent = new Event("change", {bubbles: true});
      for (const checkbox of table.getElementsByClassName("row-selector")) {
        if (checkbox.checked !== evt.target.checked) {
          checkbox.checked = evt.target.checked;
          checkbox.dispatchEvent(changeEvent);
        }
      }
    })
  } else {
    checkbox.classList.add("row-selector");
    checkbox.setAttribute("tox-family", family);
    checkbox.setAttribute("tox-gene", gene);
    checkbox.addEventListener("change", evt => {
      const tr = evt.target.closest("tr");
      tr.classList.toggle("selected", evt.target.checked);
      if (!evt.target.checked) {
        evt.target.closest("table").querySelector("#selectAll").checked = false;
      }
    })
  }
  const element = createElement(selectAll ? "th": "td", {
    children: [checkbox],
    classes: ["center", "clickable"]
  });

  element.addEventListener("click", function (evt) {
    if (evt.target === this) {
      this.firstChild.click();
    }
  });

  return element;
}

function getDetailsTableDataMap() {
  const dataMap = {};

  const tissues = dataHandler.tissues;
  const links = {
    family: {
      title: "Family",
      data(geneData, familyIdx, geneIdx) {
        return createLinkElement(geneData.family, () => createSingleDetailsTable(geneData, familyIdx, geneIdx, [
          { title: "Identifier", data() {return geneData.family} },
          {
            title: "Genes",
            data() {
              return createLinkElement(`Inspect ${dataHandler.getGeneCount(familyIdx)} members`, () => createMasterTable(
                dataHandler.genes(familyIdx).map(gene => ({family: familyIdx, gene})),
                (family, gene) => dataHandler.getGeneData(family, gene),
                dataMap.Gene
              ))
            }
          },
          { title: "Description", data() {return "..."}},
        ]));
      }
    },
    gene: {
      title: "Gene",
      data(geneData, familyIdx, geneIdx) {
        return createLinkElement(geneData.genes, () => createSingleDetailsTable(geneData, familyIdx, geneIdx, [
          { title: "Identifier", data() {return geneData.genes } },
          links.family,
          { title: "Species", data() {return geneData.species} },
          { title: "Description", data() {return "..."} },
          ...tissues.map((tissue, i) => {
            return { title: tissue, data() {return geneData.coordinates[i]} }
          })
        ]));
      }
    }
  }
  function tissueRelatedHeader(key) {
    function createValueElement(tissueData) {
      const valueElement = createElement("span");
      valueElement.innerText = tissueData[tissues.indexOf(config.get(key))].toFixed(3);
      config.onChange(key, function ({ value }) {
        valueElement.innerText = tissueData[tissues.indexOf(value)].toFixed(3);
      });
      return valueElement;
    }
    return {
      get title() {
        const title = createElement("span");
        title.innerText = config.get(key);
        config.onChange(key, ({ value }) => title.innerText = value);
        return title;
      },
      data(geneData) {
        return createValueElement(geneData.coordinates);
      }
    };
  }
  const tissueRelated = [
    tissueRelatedHeader("tissueX"),
    tissueRelatedHeader("tissueY"),
    tissueRelatedHeader("tissueZ")
  ];

  dataMap.Gene = [
    links.gene,
    links.family,
    ...tissueRelated
  ];
  dataMap.ShiftVector = [
    links.gene,
    links.family
  ];
  dataMap.Centroid = [
    links.family,
    ...tissueRelated
  ];
  return dataMap;
}

function createDataTable(options) {
  return createElement("table", { ...options, classes: ["datatable", "textselect", ...(options.classes ?? [])] });
}

function createMasterTable(elements, getData, headerMap, bodyOnly=false) {
  const tbody = createElement("tbody");

  for (const { family, gene } of elements) {
    const row = createElement("tr");
    row.appendChild(createRowSelector(false, family, gene, "Gene"));

    const data = getData(family, gene);
    for (const cell of headerMap) {
      const td = createElement("td", {
        children: [cell.data(data, family, gene)]
      });
      row.appendChild(td);
    }

    tbody.appendChild(row);
  }

  if (bodyOnly) {
    return tbody;
  } else {
    const headerCells = headerMap.map(({ title }) => createElement("th", { children: [title] }))
    const thead = createElement("thead", {
      children: [createElement("tr", { children: [createRowSelector(true), ...headerCells] })]
    });

    const table = createDataTable({ children: [thead, tbody]});
    let selectedRowCount = 0;
    table.addEventListener("change", evt => {
      if (evt.target.classList.contains("row-selector")) {
        if (evt.target.checked) {
          selectedRowCount++;
        } else {
          selectedRowCount--;
        }
        console.log(selectedRowCount)
      }
    });

    return table;
  }
}

function createSingleDetailsTable(geneData, familyIdx, geneIdx, headerMap) {
  const tBody = createElement("tbody");
  for (const header of headerMap) {
    const tr = createElement("tr", {
      children: [
        createElement("td", { children: [header.title] }),
        createElement("td", { children: [header.data(geneData, familyIdx, geneIdx)] })
      ]
    });
    tBody.appendChild(tr);
  }

  return createDataTable({ children: [tBody]});
}

function createLinkElement(value, linkContent) {
  const a = createElement("a", {
    textContent: value,
    classes: ["clickable"]
  });
  a.addEventListener("click", evt => {
    evt.preventDefault();
    show(linkContent());
  });
  return a;
}

function show(...contents) {
  const dialog = createElement("dialog", { children: contents });
  dialog.addEventListener("close", () => dialog.remove());
  document.getElementById("UI")?.appendChild(dialog);
  dialog.showModal();
}

function createUIdiv() {
  const flyer = document.getElementById("flyer");
  if (flyer) {
    const canvas = createElement("canvas", { id: "view" });
    flyer.appendChild(canvas);

    document.getElementById("UI")?.remove();
    const UI = createElement("div", { id: "UI" });
    for (const area of ["top-left", "top-right", "mid-right"]) {
      const element = createElement("aside", { id: area });
      UI.appendChild(element);
    }
    flyer.appendChild(UI);
  } else {
    throw new Error('Missing <main id="flyer"></main>');
  }
}

function createDetailButton() {
  if (!document.getElementById("detailBtn")) {
    const button = createElement("button", {
      id: "detailBtn",
      textContent: "Detail",
      classes: ["clickable"]
    });
    button.addEventListener("click", () => document.getElementById("pickedDetails")?.showModal());
    document.getElementById("mid-right")?.appendChild(button);
  }
}

function removeDetailButton() {
  document.getElementById("detailBtn")?.remove();
}

export function createTooltip(x, y, text) {
  let tooltip = document.getElementById("tooltip");
  if (tooltip === null) {
    tooltip = createElement("aside", {
      id: "tooltip",
      classes: ["tooltip"]
    });
    document.getElementById("UI")?.appendChild(tooltip);
  }
  tooltip.style.top = y + 10 + "px";
  tooltip.style.left = x + 10 + "px";
  tooltip.innerHTML = text;
}

export function removeTooltip() {
  document.getElementById("tooltip")?.remove();
}
window.addEventListener("resize", () => {
  removeTooltip();
})

function createElement(tag, options={}) {
  const element = document.createElement(tag);
  const { children, classes, ...attributes } = options;
  if (children) {
    element.append(...children);
  }
  if (classes) {
    element.classList.add(...classes);
  }
  for (const [key, value] of Object.entries(attributes)) {
    if (element[key] === undefined) {
      element.setAttribute(key, value);
    } else {
      element[key] = value;
    }
  }
  return element;
}