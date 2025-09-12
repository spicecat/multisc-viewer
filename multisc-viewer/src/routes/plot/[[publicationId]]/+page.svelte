<script lang="ts">
	import PlotForm from '$lib/components/PlotForm/PlotForm.svelte';
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
	<section class="max-w-256 mx-auto">
		<h2 class="h2 text-center">{data.publication.title}</h2>
	</section>
{/if}
<section
	use:dndzone={{ items: datasets }}
	onconsider={(e) => (datasets = e.detail.items)}
	onfinalize={(e) => (datasets = e.detail.items)}
	class="justify-center-safe flex gap-2 overflow-scroll"
>
	{#each data.plotParams.datasets as ds (`ds-${ds}`)}
		<div class="card preset-filled-surface-500 text-center">
			<div class="font-bold">{ds}</div>
			{#each data.plotParams.genes as gene (`gene-${ds}-${gene}`)}
				{#each data.plotParams.plotTypes as pt (`pt-${ds}-${pt}`)}
					{@const plotId = data.plotIds[ds][gene][pt]}
					<div class="w-xs flex aspect-video items-center justify-center p-2">
						{#await data.plots[plotId]}
							<ProgressRing value={null} />
						{:then plotData}
							<img src={plotData} alt={plotId} />
						{:catch error}
							<p>Error loading plot: {error.message}</p>
						{/await}
					</div>
				{/each}
			{/each}
		</div>
	{/each}
</section>

<hr class="hr" />
<section class="size-fit">
	<PlotForm
		datasets={data.publication?.datasets ?? data.datasets}
		genes={data.genes}
		degs={data.degs}
	/>
</section>
