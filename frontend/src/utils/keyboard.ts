/** Returns true when focus is on an input element where shortcuts should not fire. */
export function isInputFocused(target: EventTarget | null): boolean {
	if (!(target instanceof HTMLElement)) return false;
	return (
		target.tagName === "INPUT" ||
		target.tagName === "TEXTAREA" ||
		target.isContentEditable
	);
}
