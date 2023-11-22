// Taken from https://codepen.io/netsi1964/pen/JqBLPK

class ManipulateElement {
  borderTresholdMax = 20;
  borderTresholdMin = 2;

  constructor(selectorOrElement) {
    this.allowNegative = false;
    this.id = parseInt(Math.random() * 10e7 + new Date().getMilliseconds(), 10);
    this.element =
      typeof selectorOrElement === "string"
        ? document.querySelector(selectorOrElement)
        : selectorOrElement;
    if (this.element) {
      this.initCSS();
      const temp = this.element.getBoundingClientRect();
      this.dimensions = {
        left: temp.left,
        top: temp.top,
        width: temp.width,
        height: temp.height,
      };
      this.direction = null;
      this.mousePos = null;

      this.mouseMove.bind(this);
      this.element.addEventListener("pointerdown", (evt) => {
        if (evt.target.tagName.toLowerCase() === "button") {
          return;
        }
        this.pointerId = evt.pointerId;
        this.element.setPointerCapture(this.pointerId);
        this.direction = null;
      });
      this.element.addEventListener("pointerup", (evt) => {
        this.element.releasePointerCapture(this.pointerId);
        this.pointerId = null;
        this.direction = null;
        this.element.style.cursor = "grab";
      });
      this.element.addEventListener("mousemove", (evt) => {
        evt.preventDefault();
        this.mouseMove(evt);
      });
      this.endAltering.bind(this);
      this.element.addEventListener("touchend", this.endAltering);
      this.element.addEventListener("touchmove", (evt) => {
        var touch = evt.changedTouches[0];
        evt = {
          clientX: touch.clientX,
          clientY: touch.clientY,
          buttons: 1,
          shiftKey: false,
        };
        this.mouseMove(evt);
        evt.preventDefault();
        return false;
      });
    } else {
      console.log(`ManipulateElement: Could not find element "${selector}"`);
    }
  }

  initCSS() {
    var isTouchDevice =
      "ontouchstart" in window ||
      navigator.MaxTouchPoints > 0 ||
      navigator.msMaxTouchPoints > 0;
    this.element.classList.add("ManipulateElement");
    let styles = document.getElementById("ManipulateElement");

    if (styles === null) {
      const temp = document.createElement("style");
      temp.setAttribute("id", "ManipulateElement");
      temp.innerHTML = `
        ${isTouchDevice ? "* { user-select: none;}" : ""}
        body.altering .ManipulateElement {
          pointer-events: none;
        }
        body.altering .ManipulateElement.not-locked {
          pointer-events: all;
        }
              .touch-device .info-touch { ${
                isTouchDevice
                  ? "display: block; color: red; font-weight: bold;"
                  : "display: none;"
              }}
        `;
      styles = document.body.appendChild(temp);
      document.body.classList.add("touch-device");
    }

    this.useMargin =
      ["absolute", "relative", "fixed"].indexOf(
        getComputedStyle(this.element).position.toLowerCase(),
      ) === -1;
  }

  change(property, value) {
    this.element.style[property], `${value}px`;
  }

  mouseMove({ clientX, clientY, buttons, shiftKey }) {
    const BCR = this.element.getBoundingClientRect();
    const deltaX = Math.abs(clientX - BCR.x);
    const deltaY = Math.abs(clientY - BCR.y);
    const distX = BCR.width - deltaX;
    const distY = BCR.height - deltaY;

    this.direction = this.direction
      ? this.direction
      : {
          n: deltaY < this.borderTresholdMax,
          s: distY < this.borderTresholdMax,
          w: deltaX < this.borderTresholdMax,
          e: distX < this.borderTresholdMax,
        };

    const cssClass = Object.keys(this.direction).filter(
      (key) => this.direction[key],
    );

    let cursor = cssClass.join("");
    // Add the oppersite direction based on discovered direction
    switch (cursor) {
      case "e":
        cursor = "ew";
        break;
      case "w":
        cursor = "ew";
        break;
      case "n":
        cursor = "ns";
        break;
      case "s":
        cursor = "ns";
        break;
      case "ne":
        cursor = "nesw";
        break;
      case "sw":
        cursor = "nesw";
        break;
      case "nw":
        cursor = "nwse";
        break;
      case "se":
        cursor = "nwse";
        break;
      default:
        cursor = "";
        break;
    }

    this.element.style.cursor = `${
      cursor ? cursor + "-resize" : buttons === 1 ? "move" : "grab"
    }`;

    if (buttons === 1) {
      document.body.classList.add("altering");
      this.element.classList.add("not-locked");
      const changeX = this.mousePos ? this.mousePos.clientX - clientX : 0;
      const changeY = this.mousePos ? this.mousePos.clientY - clientY : 0;
      this.mousePos = { clientX, clientY };
      const closeToAnEdge = Object.keys(this.direction).filter(
        (key) => this.direction[key],
      );
      const newHeight =
        this.dimensions.height +
        (this.direction.s ? -changeY : changeY) * (shiftKey ? 2 : 1);
      const newTop = this.dimensions.top - changeY;
      const newWidth =
        this.dimensions.width +
        (this.direction.e ? -changeX : changeX) * (shiftKey ? 2 : 1);
      const newLeft = this.dimensions.left - changeX;
      const changes = [];

      const prefix = this.useMargin ? "margin-" : "";

      if (this.direction.e) {
        if (newWidth > 0) {
          if (shiftKey) {
            this.dimensions.left = newLeft + changeX * 2;
            changes.push({
              prop: `${prefix}left`,
              value: newLeft + changeX * 2,
            });
          }
          this.dimensions.width = newWidth;
          changes.push({ prop: "width", value: newWidth });
        }
      }
      if (this.direction.n) {
        if (newHeight > 0) {
          this.dimensions.height = newHeight;
          changes.push({ prop: "height", value: newHeight });
          if (this.allowNegative || newTop > 0) {
            this.dimensions.top = newTop;
            changes.push({ prop: `${prefix}top`, value: newTop });
          }
        }
      }
      if (this.direction.w) {
        if (newWidth > 0) {
          this.dimensions.width = newWidth;
          changes.push({ prop: "width", value: newWidth });

          if (this.allowNegative || newLeft > -1) {
            this.dimensions.left = newLeft;
            changes.push({ prop: `${prefix}left`, value: newLeft });
          }
        }
      }
      if (this.direction.s) {
        if (newHeight > 0) {
          if (shiftKey) {
            this.dimensions.top = newTop + changeY * 2;
            changes.push({ prop: `${prefix}top`, value: newTop + changeY * 2 });
          }
          this.dimensions.height = newHeight;
          changes.push({ prop: "height", value: newHeight });
        }
      }

      if (closeToAnEdge.length === 0) {
        // Move rectangle
        const newLeft = this.dimensions.left - changeX;
        if (this.allowNegative || newLeft > -1) {
          this.dimensions.left = newLeft;
          changes.push({ prop: `${prefix}left`, value: newLeft });
        }
        if (this.allowNegative || newTop > -1) {
          this.dimensions.top = newTop;
          changes.push({ prop: `${prefix}top`, value: newTop });
        }
      }
      const element = this;
      changes.forEach((change) => {
        this.element.style[change.prop] = `${change.value}px`;
      });
    } else {
      this.endAltering();
    }

    this.mousePos = { clientX, clientY };
  }

  endAltering() {
    document.body.classList.remove("altering");
    this.element.classList.remove("not-locked");
  }
}

export { ManipulateElement };

//   const green = new ManipulateElement(".green");
//   const blue = new ManipulateElement(".blue");
//   const img = new ManipulateElement(".img");
//   const info = new ManipulateElement(document.querySelector(".info"));
