import { publications } from '$lib/server/data';
import { type RequestHandler, json } from '@sveltejs/kit';

export const GET: RequestHandler = async () => json(Array.from(Object.values(publications)));
