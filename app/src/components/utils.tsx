export function setTheme(setColorMode: any): void {
    // check if there is a cookie
    const theme = localStorage.getItem("theme");
    if (theme) {
      document.documentElement.setAttribute("data-bs-theme", theme);
      setColorMode(theme);
    } else {
      const prefersDarkScheme = window.matchMedia("(prefers-color-scheme: dark)");
      const theme = prefersDarkScheme.matches ? "dark" : "light";
      document.documentElement.setAttribute("data-bs-theme", theme);
    setColorMode(theme);
    }
}