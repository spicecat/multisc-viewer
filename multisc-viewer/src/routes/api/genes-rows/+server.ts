import { json, type RequestHandler } from "@sveltejs/kit";
import { getGenesRows } from "$lib/server/data";

export const GET: RequestHandler = async ({ url }) => {
	const ds = url.searchParams.getAll("ds");
	const limit = Number(url.searchParams.get("limit")) || undefined;
	const offset = Number(url.searchParams.get("offset")) || undefined;
	return json(await getGenesRows(ds, limit, offset));
};
