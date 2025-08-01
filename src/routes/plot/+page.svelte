<script lang="ts">
	import ChartDisplay from '$lib/components/ChartDisplay.svelte';
	import GeneControlsDrawer from '$lib/components/GeneControls/GeneControlsDrawer.svelte';
	import GeneControlsDrawerToggle from '$lib/components/GeneControls/GeneControlsDrawerToggle.svelte';
	import html2canvas from 'html2canvas';
	import { tick } from 'svelte';
	import { dndzone } from 'svelte-dnd-action';
	import type { Gene, Grouping, Dataset } from '$lib/types/data';

	interface RenderResult {
		clustering: string;
		violin: string;
	}
	interface CompareProps {
		genes: string[];
		gene: string;
		order: string[];
	}

	const { genes, gene, order }: CompareProps = $props();
	const { query } = meta;

	let groupBy: Grouping = $state('Genotype');
	let selectedGene: string = $state(gene);
	let dsOrder = $state(order.map((ds) => ({ id: ds })));
	let bigView: boolean = $state(false);
	let bigViewCharts: (Promise<RenderResult> | null)[] = $state([null, null]);
	let geneControlsOpen = $state(false);
	let boardElement: HTMLElement | null = $state(null);
	let isDownloading = $state(false);

	const splitBy = $derived(groupBy === 'Genotype' ? 'CellType' : 'Genotype');
	const config = $derived({ selectedGene, groupBy, splitBy });

	$effect(() => {
		if (
			!history.state ||
			history.state.selectedGene !== selectedGene ||
			history.state.groupBy !== groupBy
		)
			history.pushState(
				{ selectedGene, groupBy },
				'',
				`/plot?datasets=${query.datasets}&gene=${selectedGene}&groupBy=${groupBy}&splitBy=${splitBy}`
			);
	});

	async function downloadBoard() {
		if (!boardElement || isDownloading) return;

		isDownloading = true;
		await tick();

		try {
			const canvas = await html2canvas(boardElement, {
				useCORS: true,
				allowTaint: true,
				scrollX: 0,
				scrollY: 0,
				windowWidth: boardElement.scrollWidth,
				windowHeight: boardElement.scrollHeight
			});

			const image = canvas.toDataURL('image/png');

			const link = document.createElement('a');
			link.href = image;
			link.download = `MultiSC-Viewer-${selectedGene}-${groupBy}.png`;
			document.body.appendChild(link);
			link.click();
			document.body.removeChild(link);
		} catch (error) {
			console.error('Error downloading board:', error);
		} finally {
			isDownloading = false;
		}
	}
</script>

<svelte:head>
	<title>MultiSC-Viewer - Compare</title>
</svelte:head>

<GeneControlsDrawer {genes} bind:selected={selectedGene} bind:groupBy {geneControlsOpen} />
<div class="flex-1">
	<div class="flex items-center gap-4">
		<GeneControlsDrawerToggle bind:geneControlsOpen />
		<button class="variant-filled btn" onclick={downloadBoard}>
			{#if isDownloading}
				<span>Downloading...</span>
			{:else}
				<span>Download</span>
			{/if}
		</button>
	</div>
	<div class="board" bind:this={boardElement} use:dndzone={{ items: dsOrder, flipDurationMs: 300 }}>
		{#each dsOrder as { id: ds }, i (ds)}
			<ChartDisplay {dataset} {config} {bigView} {bigViewCharts} />
		{/each}
	</div>
</div>

<style lang="scss">
	.board {
		min-height: 40vh;
		margin: 0.5em;
		display: flex;
		overflow-x: scroll;
		background-color: white;
		padding: 1em;
	}
</style>
