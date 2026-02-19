import { error } from "@sveltejs/kit";
import { getDatasets } from "$lib/server/data";
import type { PageServerLoad } from "./$types";

export const load: PageServerLoad = async ({ params }) => {
	const datasets = await getDatasets([params.dsId]);
	const dataset = datasets[params.dsId];
	if (!dataset) error(404, `Dataset ${params.dsId} not found`);
	const name = dataset.displayName ?? dataset.name ?? dataset.id;
	return {
		dataset,
		meta: {
			title: `MultiSC-Viewer - Dataset ${name}`,
			description: `Explore dataset "${name}" with MultiSC-Viewer, a web application for visualizing and comparing gene expression across multiple datasets, brain regions, disease conditions, and species.`,
		},
	};
};
