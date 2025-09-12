<script lang="ts">
	import { enhance } from '$app/forms';
	import { page } from '$app/state';
	import type { Dataset, DEGs, Gene } from '$lib/types/data';
	import DatasetsTable from '../DataTable/DatasetsTable.svelte';
	import GenesTable from '../DataTable/GenesTable.svelte';
	import Groupings from './Groupings.svelte';
	import PlotTypes from './PlotTypes.svelte';

	let { datasets, genes, degs }: { datasets: Dataset[]; genes?: Gene[]; degs?: DEGs } = $props();

	let publicationId = $derived(page.params.publicationId);
	let selectedDEGs = $state('');
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
	{#if degs}
		<div class="flex gap-4">
			<label class="label">
				<span class="label-text">DEGs</span>
				<select class="select" bind:value={selectedDEGs} name="degsId">
					<option value=""></option>
					{#each Object.keys(degs) as degsId (`degs-${degsId}`)}
						<option value={degsId}>{degsId}</option>
					{/each}
				</select>
				<textarea
					class="textarea"
					placeholder="Enter genes, one per line"
					value={degs[selectedDEGs]?.join('\n')}
				></textarea>
			</label>
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
