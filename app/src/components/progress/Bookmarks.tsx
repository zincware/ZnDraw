import { Tooltip } from "@mui/material";
import { FaRegBookmark } from "react-icons/fa";

export const Bookmarks = ({
	length,
	bookmarks,
	setBookmarks,
	setStep,
}: {
	length: number;
	bookmarks: { [key: number]: string };
	setBookmarks: (bookmarks: { [key: number]: string }) => void;
	setStep: (step: number) => void;
}) => {
	const handleBookmarkClick = (event: any, number: number) => {
		if (event.shiftKey) {
			const newBookmarks = { ...bookmarks };
			delete newBookmarks[number];
			setBookmarks(newBookmarks);
		} else {
			setStep(number);
		}
	};

	return (
		<>
			{Object.keys(bookmarks).map((key) => {
				const position = Number.parseInt(key);
				return (
					<Tooltip
						key={position}
						title={bookmarks[position]}
						slotProps={{
							tooltip: {
								style: { marginLeft: "0.375em" },
							},
						}}
					>
						<div
							className="position-absolute progress-bar-bookmark"
							style={{
								left: `${(position / length) * 100}%`,
								top: 7,
								marginLeft: "-0.375em",
							}}
						>
							<FaRegBookmark
								className="position-absolute"
								size={"0.75em"}
								onClick={(e) => handleBookmarkClick(e, position)}
							/>
						</div>
					</Tooltip>
				);
			})}
		</>
	);
};