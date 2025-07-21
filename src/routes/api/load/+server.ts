// POST /api/load
import { json } from '@sveltejs/kit';
import { daemonService } from '$lib/server/daemon.service';

export async function POST({ request }) {
  const { datasets } = await request.json();
  await daemonService.load(datasets);
  return json({ status: 'ok' });
}
