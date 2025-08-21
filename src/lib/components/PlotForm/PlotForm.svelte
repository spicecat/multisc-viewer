<script lang="ts">
	import { enhance } from '$app/forms';
	import { page } from '$app/state';
	import type { Dataset, Gene } from '$lib/types/data';
	import DatasetsTable from '../DataTable/DatasetsTable.svelte';
	import GenesTable from '../DataTable/GenesTable.svelte';
	import Groupings from './Groupings.svelte';
	import PlotTypes from './PlotTypes.svelte';

	let { datasets, genes }: { datasets: Dataset[]; genes?: Gene[] } = $props();

	let publicationId = $derived(page.params.publicationId);
</script>

<form
	method="POST"
	action="/plot?/plot"
	use:enhance={() =>
		async ({ update }) =>
			update({ invalidateAll: true })}
	class="space-y-4"
>
	<input type="hidden" name="publicationId" value={publicationId} />
	{#if genes}
		<div class="flex gap-4">
			<Groupings />
			<PlotTypes />
		</div>
	{/if}
	<div class="flex gap-4">
		{#if genes}
			<div class="space-y-2">
				<GenesTable {genes} />
			</div>
		{/if}
		<DatasetsTable {datasets} />
	</div>
</form>
