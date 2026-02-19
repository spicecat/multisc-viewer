import { json, type RequestHandler } from "@sveltejs/kit";
import { getDatasets } from "$lib/server/data";

export const GET: RequestHandler = async ({ url }) => {
	const ds = url.searchParams.getAll("ds");
	return json(await getDatasets(ds.length ? ds : undefined));
};
