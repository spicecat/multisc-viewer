<script lang="ts">
	import type { Dataset } from '$lib/types/data';
	import { datasetColumns } from '$lib/types/data-table';
	import { Tabs } from '@skeletonlabs/skeleton-svelte';
	import DataTableSearch from './DataTableSearch.svelte';

	let { datasets }: { datasets: Dataset[] } = $props();

	let cellType = $state('');
	const cellTypes = $derived([...new Set(datasets.map((ds) => ds.cellType))]);
</script>

<Tabs value={cellType} onValueChange={(e) => (cellType = e.value)} fluid>
	{#snippet list()}
		<Tabs.Control value="">All</Tabs.Control>
		{#each cellTypes as ct}
			<Tabs.Control value={ct}>{ct}</Tabs.Control>
		{/each}
	{/snippet}
	{#snippet content()}
		<DataTableSearch
			label="Datasets"
			data={cellType ? datasets.filter((ds) => ds.cellType === cellType) : datasets}
			columns={datasetColumns}
			select="checkbox"
		/>
	{/snippet}
</Tabs>
