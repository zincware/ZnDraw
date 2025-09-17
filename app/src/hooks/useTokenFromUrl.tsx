import { useEffect, useState } from "react";
import { useSearchParams } from "react-router-dom";

export const useTokenFromUrl = () => {
	const [searchParams] = useSearchParams();
	const [token, setToken] = useState<string | null>(null);

	useEffect(() => {
		const urlToken = searchParams.get("token");
		if (urlToken) {
			setToken(urlToken);
			console.log("Token extracted from URL:", urlToken);
		}
	}, [searchParams]);

	return token;
};
