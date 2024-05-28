class Bookmarks {
  constructor(world, cache, socket) {
    this.world = world;
    this.cache = cache;
    this.socket = socket;
    this.bookmarks = {};

    // on keypress b set the bookmark
    window.addEventListener("keypress", (event) => {
      if (document.activeElement === document.body && event.key === "b") {
        // get the current step
        const step = this.world.getStep();
        // // set the bookmark at the current step
        this.bookmarks[step] = `Bookmark ${step}`;
        // // update the bookmarks
        this.updateBookmarks();
        this.socket.emit("room:bookmarks:set", this.bookmarks);
      }
    });
  }

  updateBookmarks() {
    const bookmark_envelope = document.getElementById("bookmarks");
    // remove all children
    while (bookmark_envelope.firstChild) {
      bookmark_envelope.removeChild(bookmark_envelope.firstChild);
    }

    // remove all bookmarks for frames larger than the cache length
    for (const [index, name] of Object.entries(this.bookmarks)) {
      if (index > this.cache.get_length()) {
        delete this.bookmarks[index];
      }
    }

    for (const [index, name] of Object.entries(this.bookmarks)) {
      const button = document.createElement("button");
      button.type = "button";
      button.className = "btn btn-link";
      button.innerHTML = '<i class="fa-regular fa-bookmark"></i>';
      // add data-bs-toggle="tooltip" data-bs-placement="top" title=name
      button.setAttribute("data-bs-placement", "top");
      button.setAttribute("data-bs-title", name);
      // set z-index to 100
      button.style.zIndex = "100";
      button.addEventListener("click", (event) => {
        // check if shift key is being pressed
        if (event.shiftKey) {
          // delete the bookmark
          delete this.bookmarks[index];
          this.updateBookmarks();
        } else {
          this.world.setStep(index);
        }
      });
      // set the position of the button to be relative  index / this.cache.get_length()
      button.style.left = `${(index / this.cache.get_length()) * 100}%`;

      button.style.position = "absolute";
      button.style.bottom = "5px";
      bookmark_envelope.appendChild(button);

      // set the tooltip
      const tooltip = new bootstrap.Tooltip(button, {
        boundary: document.body,
      });
      // hide the tooltip on click
      button.addEventListener("click", (event) => {
        tooltip.hide();
      });
    }
  }

  set(bookmarks) {
    this.bookmarks = bookmarks;
    this.updateBookmarks();
  }

  step() {
    this.updateBookmarks();
    return true;
  }

  jumpNext() {
    const step = this.world.getStep();
    const bookmarks = Object.keys(this.bookmarks);
    const next_bookmark = bookmarks.find((bookmark) => bookmark > step);
    if (next_bookmark) {
      this.world.setStep(next_bookmark);
    }
  }

  jumpPrevious() {
    const step = this.world.getStep();
    const bookmarks = Object.keys(this.bookmarks);
    const previous_bookmark = bookmarks
      .reverse()
      .find((bookmark) => bookmark < step);
    if (previous_bookmark) {
      this.world.setStep(previous_bookmark);
    }
  }
}

export { Bookmarks };
