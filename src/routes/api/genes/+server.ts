// GET /api/genes?datasets=...
import { json } from '@sveltejs/kit';
import { dataService } from '$lib/server/data.service';

export async function GET({ url }) {
  const datasets = url.searchParams.get('datasets')?.split(',') || [];
  return json(dataService.getGenes(datasets));
}
