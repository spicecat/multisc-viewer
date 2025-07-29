<script lang="ts">
	import DataTableSearch from './DataTableSearch.svelte';

	const { publications }: { publications: Publication[] } = $props();

	const publicationData = $derived(
		publications.map((publication) => ({
			...publication,
			publication: {
				name: publication.title,
				href: `/publication/${publication.publicationId}`
			},
			datasets: publication.datasets,
			// year: [...new Set(publication.datasets.map((ds) => ds.year))],
			// PMID: [...new Set(publication.datasets.map((ds) => ds.PMID))],
			species: [...new Set(publication.datasets.map((ds) => ds.species))],
			// author: [...new Set(publication.datasets.map((ds) => ds.author))],
			disease: [...new Set(publication.datasets.map((ds) => ds.disease).flat())],
			cellType: [...new Set(publication.datasets.map((ds) => ds.cellType))]
		}))
	);

	const publicationColumns: Column[] = [
		{ key: 'publication', label: 'Publication', url: true },
		{ key: 'author', label: 'Authors' },
		{ key: 'journal', label: 'Journal' },
		{ key: 'year', label: 'Year' },
		{ key: 'PMID', label: 'PMID' },
		{ key: 'disease', label: 'Disease' },
		{ key: 'cellType', label: 'Cell Type' },
		{ key: 'datasets', label: 'Datasets' }
	];
</script>

<DataTableSearch label="Publications" data={publicationData} columns={publicationColumns} />
