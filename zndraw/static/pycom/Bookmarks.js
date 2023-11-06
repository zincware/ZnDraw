class Bookmarks {
  constructor(world, cache, socket) {
    this.world = world;
    this.cache = cache;
    this.socket = socket;
    this.bookmarks = {};

    this.socket.on(
      "bookmarks:get",
      function (callback) {
        console.log("bookmarks:get");
        callback(this.bookmarks);
      }.bind(this),
    );
    this.socket.on("bookmarks:set", (bookmarks) => {
      console.log("bookmarks:set");
      this.bookmarks = bookmarks;
      this.updateBookmarks();
    });

    // on keypress b set the bookmark
    window.addEventListener("keypress", (event) => {
      if (document.activeElement === document.body && event.key === "b") {
        // get the current step
        const step = this.world.getStep();
        // // set the bookmark at the current step
        this.bookmarks[step] = `Bookmark ${step}`;
        // // update the bookmarks
        this.updateBookmarks();
      }
    });
  }

  updateBookmarks() {
    const bookmark_envelope = document.getElementById("bookmarks");
    // remove all children
    while (bookmark_envelope.firstChild) {
      bookmark_envelope.removeChild(bookmark_envelope.firstChild);
    }
    console.log("update bookmarks");

    for (const [index, name] of Object.entries(this.bookmarks)) {
      const button = document.createElement("button");
      button.type = "button";
      button.className = "btn btn-link";
      button.innerHTML = '<i class="fa-regular fa-bookmark"></i>';
      // add data-bs-toggle="tooltip" data-bs-placement="top" title=name
      button.setAttribute("data-bs-toggle", "tooltip");
      button.setAttribute("data-bs-placement", "top");
      button.setAttribute("title", name);
      // set z-index to 100
      button.style.zIndex = "100";
      button.addEventListener("click", () => {
        this.world.setStep(index);
      });
      // set the position of the button to be relative  index / this.cache.get_length()
      button.style.left = `${(index / this.cache.get_length()) * 100}%`;

      button.style.position = "absolute";
      button.style.bottom = "5px";
      console.log(button.style.left);
      bookmark_envelope.appendChild(button);
    }
  }

  step() {
    this.updateBookmarks();
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
