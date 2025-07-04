"use strict";

export function createTooltip(x, y, text) {
  let tooltip = document.getElementById("tooltip");
  if (tooltip === null) {
    tooltip = document.createElement("aside");
    tooltip.id = "tooltip";
    tooltip.classList.add("tooltip");
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