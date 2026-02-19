import { json, type RequestHandler } from "@sveltejs/kit";
import { getDEGs } from "$lib/server/data";

export const GET: RequestHandler = async ({ url }) => {
	const datasets = url.searchParams.getAll("ds");
	return json(await getDEGs(datasets));
};
