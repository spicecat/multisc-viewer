<script lang="ts">
import { Camera } from "@lucide/svelte";
import { Progress } from "@skeletonlabs/skeleton-svelte";
import html2canvas from "html2canvas-pro";
import { SvelteMap } from "svelte/reactivity";
import { dndzone } from "svelte-dnd-action";
import Publication from "$lib/components/Data/Publication.svelte";
import PlotsForm from "$lib/components/Plots/Form.svelte";
import type { PageProps } from "./$types";

const { data }: PageProps = $props();
const {
	datasets,
	publication,
	genesRows,
	plotsIdMap,
	plotsParams,
	plotsResults,
} = $derived(data);

const plotsElements = new SvelteMap<string, HTMLImageElement>();
async function downloadPlots() {
	const plotIds = Object.keys(plotsResults);
	for (const [plotId, element] of plotsElements) {
		if (!plotIds.includes(plotId) || !element) continue;
		try {
			const canvas = await html2canvas(element);
			const image = canvas.toDataURL("image/png");
			const link = document.createElement("a");
			link.download = `${plotId}.png`;
			link.href = image;
			link.click();
			link.remove();
		} catch (error) {
			console.error("Error capturing screenshot:", error);
		}
	}
}

let items = $derived(plotsParams.ds.map((d) => ({ id: d })));
</script>

<svelte:head>
	<title>{data.meta.title}</title>
	<meta name="description" content={data.meta.description} />
</svelte:head>

{#if publication}
	<Publication {publication} />
	<hr class="hr" />
{/if}

<section class="space-y-2">
	<button type="button" class="btn preset-filled" onclick={downloadPlots}>
		Download Plots <Camera size={18} />
	</button>

	<section
		use:dndzone={{ items }}
		onconsider={(e) => (items = e.detail.items)}
		onfinalize={(e) => (items = e.detail.items)}
		class="flex justify-center-safe gap-2 overflow-scroll"
	>
		{#each items as ds (ds.id)}
			{@const d = ds.id}
			<div class="card preset-filled-surface-500 text-center">
				<div class="font-bold">{datasets[d]?.displayName}</div>
				{#each plotsParams.gene as g (`plots-gene-${d}-${g}`)}
					{#each plotsParams.pt as p (`plots-pt-${d}-${g}-${p}`)}
						{@const plotId = plotsIdMap[d]?.[g]?.[p]}
						{#if plotId}
							<div class="flex w-xs flex-col items-center p-2">
								{#await plotsResults[plotId]}
									<Progress value={null} class="items-center">
										<Progress.Circle>
											<Progress.CircleTrack />
											<Progress.CircleRange />
										</Progress.Circle>
									</Progress>
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
</section>

<hr class="hr" />

<section class="mx-auto">
	<PlotsForm {datasets} {genesRows} />
</section>
