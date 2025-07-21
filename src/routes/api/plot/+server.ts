// GET /api/plot?dataset=...&gene=...&groupBy=...&splitBy=...
import { json } from '@sveltejs/kit';
import { plotService } from '$lib/server/plot.service';

export async function GET({ url }) {
  const dataset = url.searchParams.get('dataset');
  const gene = url.searchParams.get('gene');
  const groupBy = url.searchParams.get('groupBy');
  const splitBy = url.searchParams.get('splitBy');
  const result = await plotService.render(dataset, gene, groupBy, splitBy);
  return json(result);
}
