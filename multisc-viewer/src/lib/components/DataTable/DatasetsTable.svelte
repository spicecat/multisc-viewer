<script lang="ts">
	import { page } from '$app/state';
	import type { Dataset } from '$lib/types/data';
	import { type DatasetData, datasetColumns } from '$lib/types/data-table';
	import { Tabs } from '@skeletonlabs/skeleton-svelte';
	import { omit, uniq } from 'lodash-es';
	import DataTableSearch from './DataTableSearch.svelte';

	let { datasets }: { datasets: Dataset[] } = $props();

	let cellType = $state('');
	const cellTypes = $derived(uniq(datasets.map((ds) => ds.cellType)));

	const data: DatasetData[] = $derived(
		(cellType ? datasets.filter((ds) => ds.cellType === cellType) : datasets)
			.map((ds) => omit(ds, ['defaultGenes', 'size']))
			.map((ds) => ({
				...ds,
				pubmed: new URL(`https://pubmed.ncbi.nlm.nih.gov/${ds.PMID}`)
			}))
	);
	const select = $derived(page.url.searchParams.getAll('ds'));
</script>

<DataTableSearch name="datasets" {data} columns={datasetColumns} {select}>
	<div class="flex items-center gap-2">
		<button formaction="/plot?/plot" type="submit" class="btn preset-filled">Plot</button>
		<Tabs value={cellType} onValueChange={(e) => (cellType = e.value)} fluid>
			{#snippet list()}
				<Tabs.Control value="">All</Tabs.Control>
				{#each cellTypes as ct (`tab-${ct}`)}
					<Tabs.Control value={ct}>{ct}</Tabs.Control>
				{/each}
			{/snippet}
		</Tabs>
	</div>
</DataTableSearch>
