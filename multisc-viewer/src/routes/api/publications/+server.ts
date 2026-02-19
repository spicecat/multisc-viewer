import { json, type RequestHandler } from "@sveltejs/kit";
import { getPublications } from "$lib/server/data";

export const GET: RequestHandler = async ({ url }) => {
	const pub = url.searchParams.getAll("pub");
	return json(await getPublications(pub));
};
