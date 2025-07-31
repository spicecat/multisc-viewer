<script lang="ts">
	import DatasetsTable from '$lib/components/DataTable/DatasetsTable.svelte';
	import type { Dataset } from '$lib/types/data';
	import Button, { Label } from '@smui/button';

	let { datasets }: { datasets: Dataset[] } = $props();

	let selectedDatasets: string[] = $state([]);
	let href = $derived(
		`/plot?${new URLSearchParams(selectedDatasets.map((ds) => ['ds', ds])).toString()}`
	);
</script>

<div style="display: flex;">
	<div style="align-items: center;display: flex;gap: 1rem;">
		<Button variant="raised" {href} {...{ disabled: Boolean(selectedDatasets.length) }}>
			<Label>Plot</Label>
		</Button>
	</div>

	<DatasetsTable {datasets} bind:selected={selectedDatasets} />
</div>
