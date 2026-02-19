import { getDatasets, getPublications } from "$lib/server/data";
import type { PageServerLoad } from "./$types";

export const load: PageServerLoad = async () => {
	const publications = await getPublications();
	const ds = Object.values(publications).flatMap((pub) => pub.datasets);
	const datasets = await getDatasets(ds);
	return {
		datasets,
		publications,
		meta: {
			title: "MultiSC-Viewer - Publications",
			description:
				"Explore publications with MultiSC-Viewer, a web application for visualizing and comparing gene expression across multiple datasets, brain regions, disease conditions, and species.",
		},
	};
};
