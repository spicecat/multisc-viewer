<script lang="ts">
	import type { Dataset } from '$lib/types/data';
	import Tab, { Label } from '@smui/tab';
	import TabBar from '@smui/tab-bar';
	import DataTableSearch from './DataTableSearch.svelte';
	import { datasetColumns } from '$lib/types/data-table';

	let {
		datasets,
		selected = $bindable()
	}: {
		datasets: Dataset[];
		selected: string[];
	} = $props();

	let cellType = $state('');

	let cellTypes = $derived([...new Set(datasets.map((ds) => ds.cellType))]);
	let cellDatasets = $derived(datasets.filter((ds) => ds.cellType === cellType));

	$effect(() => {
		if (cellTypes.length > 0 && !cellType) cellType = cellTypes[0];
	});
</script>

<TabBar tabs={cellTypes} bind:active={cellType}>
	{#snippet tab(tab)}
		<Tab>
			<Label>{tab}</Label>
		</Tab>
	{/snippet}
</TabBar>

<DataTableSearch label="Datasets" data={cellDatasets} columns={datasetColumns} bind:selected />
