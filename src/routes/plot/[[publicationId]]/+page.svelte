<script lang="ts">
	import { Progress } from '@skeletonlabs/skeleton-svelte';
	import { dndzone } from 'svelte-dnd-action';
	import type { PageData } from './$types';

	let { data }: { data: PageData } = $props();
	let datasets = $state(data.datasets);
</script>

<svelte:head>
	<title>Plot {data.datasets.map((ds) => ds.id).join(' ')}</title>
</svelte:head>

{#await data.publication}
	<p>Loading publication...</p>
{:then publication}
	<section class="mx-auto max-w-256">
		<h2 class="text-center h2">{publication.title}</h2>
	</section>
{/await}
<section
	use:dndzone={{ items: datasets }}
	onconsider={(e) => (datasets = e.detail.items)}
	onfinalize={(e) => (datasets = e.detail.items)}
	class="flex gap-2 overflow-x-scroll"
>
	{#each datasets as ds (ds.id)}
		<div class="card preset-filled-surface-500 text-center">
			<div class="preset-typo-caption">{ds.id}</div>
			{#each data.plotParams.plotTypes as pt}
				{@const plotId = data.plotIds[ds.id]?.[pt]}
				<div class="flex aspect-video w-xs items-center p-2">
					{#await data.plots[plotId]}
						<Progress value={null} />
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
