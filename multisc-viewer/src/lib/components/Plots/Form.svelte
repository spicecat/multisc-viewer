<script lang="ts">
	import { navigating, page } from '$app/state';
	import type { Datasets, GenesRows } from '$lib/types/daemon';
	import { ArrowRight } from '@lucide/svelte';
	import { Progress } from '@skeletonlabs/skeleton-svelte';
	import DatasetsTable from '../Data/DatasetsTable.svelte';
	import GenesTable from '../Data/GenesTable.svelte';

	import Groupings from './Groupings.svelte';
	import PlotTypes from './PlotTypes.svelte';

	let { datasets, genesRows }: { datasets: Datasets; genesRows?: GenesRows } = $props();
</script>

<form action="/plots" class="space-y-4">
	<input type="hidden" name="pub" value={page.params.pubId} />
	{#if genesRows}
	<div class="mx-auto flex gap-4">
		<div>
			<GenesTable {datasets} {genesRows} />
		</div>
		<div class="space-y-4">
			<Groupings />
			<PlotTypes />
		</div>
	</div>
	{/if}
	<DatasetsTable {datasets}>
		<button type="submit" class="my-auto btn flex items-center preset-filled-primary-500">
			Plot Datasets
			{#if navigating.to}
				<Progress value={null}>
					<Progress.Circle style="--size: 24px; --thickness: 4px;">
						<Progress.CircleTrack />
						<Progress.CircleRange />
					</Progress.Circle>
				</Progress>
			{:else}
				<ArrowRight size="24" />
			{/if}
		</button>
	</DatasetsTable>
</form>
