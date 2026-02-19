import { json, type RequestHandler } from "@sveltejs/kit";
import { getGenes } from "$lib/server/data";

export const GET: RequestHandler = async ({ url }) => {
	const ds = url.searchParams.getAll("ds");
	return json(await getGenes(ds));
};
