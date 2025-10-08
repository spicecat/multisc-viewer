<script lang="ts">
	import { navigating, page } from '$app/state';
	import type { Datasets, DEGs, Genes } from '$lib/types/daemon';
	import { ArrowRight } from '@lucide/svelte';
	import { ProgressRing } from '@skeletonlabs/skeleton-svelte';
	import DatasetsTable from '../Data/DatasetsTable.svelte';
	import GenesTable from '../Data/GenesTable.svelte';

	import Groupings from './Groupings.svelte';
	import PlotTypes from './PlotTypes.svelte';

	let { datasets, genes, degs }: { datasets: Datasets; genes?: Genes; degs?: DEGs } = $props();
</script>

<form action="/plots" class="space-y-4">
	<input type="hidden" name="pub" value={page.params.pubId} />
	<div class="mx-auto flex gap-4">
		{#if genes && degs}
			<div>
				<GenesTable {datasets} {genes} {degs} />
			</div>
			<div class="space-y-4">
				<Groupings />
				<PlotTypes />
			</div>
		{/if}
	</div>
	<DatasetsTable {datasets}>
		<button type="submit" class="btn preset-filled-primary-500">
			Plot Datasets
			{#if navigating.to}
				<ProgressRing value={null} size="size-5" />
			{:else}
				<ArrowRight />
			{/if}
		</button>
	</DatasetsTable>
</form>
