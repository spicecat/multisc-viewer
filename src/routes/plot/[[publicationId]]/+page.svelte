<script lang="ts">
	import DatasetsTable from '$lib/components/DataTable/DatasetsTable.svelte';
	import PlotControls from '$lib/components/PlotControls/PlotControls.svelte';
	import { ProgressRing } from '@skeletonlabs/skeleton-svelte';
	import { dndzone } from 'svelte-dnd-action';
	import type { PageProps } from './$types';

	let { data }: PageProps = $props();
	let datasets = $state(data.datasets);
</script>

<svelte:head>
	<title>Plot {data.datasets.map((ds) => ds.id).join(' ')}</title>
</svelte:head>

{#if data.publication}
	<section class="mx-auto max-w-256">
		<h2 class="text-center h2">{data.publication.title}</h2>
	</section>
{/if}

<section
	use:dndzone={{ items: datasets }}
	onconsider={(e) => (datasets = e.detail.items)}
	onfinalize={(e) => (datasets = e.detail.items)}
	class="flex justify-center gap-2 overflow-x-scroll"
>
	{#each datasets as ds (`ds-${ds.id}`)}
		<div class="card preset-filled-surface-500 text-center">
			<div class="font-bold">{ds.id}</div>
			{#each data.plotParams.plotTypes as pt (`pt-${ds.id}-${pt}`)}
				{@const plotId = data.plotIds[ds.id]?.[pt]}
				<div class="flex aspect-video w-xs items-center justify-center p-2">
					{#await data.plots[plotId]}
						<ProgressRing value={null} />
					{:then plotData}
						<img src={plotData} alt={plotId} />
					{:catch error}
						<p>Error loading plot: {error.message}</p>
					{/await}
				</div>
			{/each}
		</div>
	{/each}
</section>

<form method="POST" action="?/plot">
	<hr class="hr" />
	<section class="size-fit">
		<h2 class="text-center h2">Plot</h2>
		<div class="size-fit">
			<button type="submit" class="btn preset-filled">Plot</button>
		</div>
		<div class="flex gap-4">
			<PlotControls genes={data.genes} />
			<DatasetsTable datasets={data.publication?.datasets ?? data.datasets} />
		</div>
	</section>
</form>
