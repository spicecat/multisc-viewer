<script lang="ts">
	import Publication from '$lib/components/Data/Publication.svelte';
	import PlotsForm from '$lib/components/Plots/Form.svelte';
	import { Camera } from '@lucide/svelte';
	import { ProgressRing } from '@skeletonlabs/skeleton-svelte';
	import html2canvas from 'html2canvas-pro';
	import { dndzone } from 'svelte-dnd-action';
	import { SvelteMap } from 'svelte/reactivity';
	import type { PageProps } from './$types';

	let { data }: PageProps = $props();
	let { datasets, publication, genes, degs, plotsIdMap, plotsParams, plotsResults } =
		$derived(data);

	let plotsElements = new SvelteMap<string, HTMLImageElement>();
	async function downloadPlots() {
		const plotIds = Object.keys(plotsResults);
		for (const [plotId, element] of plotsElements) {
			if (!plotIds.includes(plotId) || !element) continue;
			try {
				const canvas = await html2canvas(element);
				const image = canvas.toDataURL('image/png');
				const link = document.createElement('a');
				link.download = `${plotId}.png`;
				link.href = image;
				link.click();
				link.remove();
			} catch (error) {
				console.error('Error capturing screenshot:', error);
			}
		}
	}

	let items = $derived(plotsParams.ds.map((d) => ({ id: d, _id: d })));
</script>

<svelte:head>
	<title
		>Plots {publication &&
			Object.values(datasets)
				.map((ds) => ds.displayName ?? ds.name ?? ds._id)
				.join(' ')}
	</title>
</svelte:head>

{#if publication}
	<Publication {publication} />
	<hr class="hr" />
{/if}

<button type="button" class="btn preset-filled" onclick={downloadPlots}>
	Download All Plots <Camera size={18} />
</button>

<section
	use:dndzone={{ items }}
	onconsider={(e) => (items = e.detail.items)}
	onfinalize={(e) => (items = e.detail.items)}
	class="flex justify-center-safe gap-2 overflow-scroll"
>
	{#each items as ds (ds.id)}
		{@const d = ds._id}
		<div class="card preset-filled-surface-500 text-center">
			<div class="font-bold">{datasets[d]?.displayName}</div>
			{#each plotsParams.gene as g (`plots-gene-${d}-${g}`)}
				{#each plotsParams.pt as p (`plots-pt-${d}-${g}-${p}`)}
					{@const plotId = plotsIdMap[d]?.[g]?.[p]}
					{#if plotId}
						<div class="flex aspect-video w-xs items-center justify-center p-2">
							{#await plotsResults[plotId]}
								<ProgressRing value={null} />
							{:then plotData}
								<img
									bind:this={() => plotsElements.get(plotId), (e) => plotsElements.set(plotId, e)}
									src={plotData}
									alt={`Plot ${plotId}`}
								/>
							{:catch error}
								<p>Error loading plot: {error.message}</p>
							{/await}
						</div>
					{/if}
				{/each}
			{/each}
		</div>
	{/each}
</section>

<hr class="hr" />

<section class="mx-auto">
	<PlotsForm {datasets} {genes} {degs} />
</section>
