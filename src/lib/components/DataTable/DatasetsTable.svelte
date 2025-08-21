<script lang="ts">
	import { enhance } from '$app/forms';
	import type { Dataset } from '$lib/types/data';
	import { type DatasetData, datasetColumns } from '$lib/types/data-table';
	import { Tabs } from '@skeletonlabs/skeleton-svelte';
	import { omit, uniq } from 'lodash-es';
	import DataTableSearch from './DataTableSearch.svelte';

	let { datasets, publicationId }: { datasets: Dataset[]; publicationId?: string } = $props();

	let cellType = $state('');
	const cellTypes = $derived(uniq(datasets.map((ds) => ds.cellType)));

	const data: DatasetData[] = $derived(
		(cellType ? datasets.filter((ds) => ds.cellType === cellType) : datasets).map((ds) =>
			omit(ds, ['size', 'defaultGene'])
		)
	);
</script>

<section class="mx-auto size-fit">
	<h2 class="text-center h2">Datasets</h2>
	<form method="POST" action="/plot?/plot" use:enhance>
		<input type="hidden" name="publicationId" value={publicationId} />
		<div class="mb-4">
			<button type="submit" class="btn preset-filled">Plot</button>
		</div>
		<DataTableSearch name="datasets" {data} columns={datasetColumns} select="checkbox">
			<Tabs value={cellType} onValueChange={(e) => (cellType = e.value)} fluid>
				{#snippet list()}
					<Tabs.Control value="">All</Tabs.Control>
					{#each cellTypes as ct (`tab-${ct}`)}
						<Tabs.Control value={ct}>{ct}</Tabs.Control>
					{/each}
				{/snippet}
			</Tabs>
		</DataTableSearch>
	</form>
</section>
