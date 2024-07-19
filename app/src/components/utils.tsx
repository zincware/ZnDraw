import { useEffect, useState } from "react";

export const useColorMode = (): [string, () => void] => {
  const [colorMode, setColorMode] = useState<string>("light");

  useEffect(() => {
    const theme =
      localStorage.getItem("theme") ||
      (window.matchMedia("(prefers-color-scheme: dark)").matches
        ? "dark"
        : "light");
    setTheme(theme, setColorMode);
  }, []);

  const handleColorMode = () => {
    const newColorMode = colorMode === "light" ? "dark" : "light";
    setTheme(newColorMode, setColorMode);
    localStorage.setItem("theme", newColorMode);
  };

  return [colorMode, handleColorMode];
};

const setTheme = (
  theme: string,
  setColorMode: (mode: string) => void,
): void => {
  document.documentElement.setAttribute("data-bs-theme", theme);
  setColorMode(theme);
};
