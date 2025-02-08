<script lang="ts">
	import ChartDisplay from '$lib/components/ChartDisplay.svelte';
	import meta from '$meta';
	import { dndzone } from 'svelte-dnd-action';
	import { self } from 'svelte/legacy';

	const { genes, gene, order }: CompareProps = $props();

	const { query } = meta;

	let geneSearch: string = $state(''),
		groupBy: string = $state('Genotype'),
		selectedGene: string = $state(gene),
		dsOrder = $state(order.map((ds) => ({ id: ds }))),
		bigView: boolean = $state(false),
		bigViewCharts: (Promise<RenderResult> | null)[] = $state([null, null]);

	const filteredGenes = $derived(genes.filter((gene) => gene.toLowerCase().includes(geneSearch.toLowerCase())).slice(0, 100)),
		splitBy = $derived(groupBy === 'Genotype' ? 'CellType' : 'Genotype'),
		config = $derived({ selectedGene, groupBy, splitBy });

	$effect(() => {
		if (!history.state || history.state.selectedGene !== selectedGene || history.state.groupBy !== groupBy) {
			history.pushState({ selectedGene, groupBy }, '', `/compare?datasets=${query.datasets}&gene=${selectedGene}&groupBy=${groupBy}&splitBy=${splitBy}`);
		}
	});
</script>

<div class="row controls">
	<div class="col center">
		<div class="nav">
			<a href="/" role="button">Go Back</a>
		</div>
		<div class="group-split">
			<div class="row">
				<div class="col">
					<h4>Group By</h4>
					<p>{groupBy}</p>
				</div>
				<div class="col swap">
					<button onclick={() => (groupBy = splitBy)}>
						<i class="fa-solid fa-left-right"></i>
					</button>
				</div>
				<div class="col">
					<h4>Split By</h4>
					<p>{splitBy}</p>
				</div>
			</div>
		</div>
		<div class="row apart">
			{#if bigView}
				<h6>Click datasets, or</h6>
				<button
					class="error"
					onclick={() => {
						bigView = false;
						bigViewCharts = [null, null];
					}}>Cancel</button
				>
			{:else}
				<button onclick={() => (bigView = true)}>Big View</button>
			{/if}
		</div>
	</div>
	<div class="gene-select">
		<label>
			<h2>Gene</h2>
			<input type="text" placeholder="Search" bind:value={geneSearch} />
		</label>
		<ul class="gene-list">
			{#each filteredGenes as gene}
				<li class="gene-option" class:selected={selectedGene === gene} onclick={() => (selectedGene = gene)}>
					{gene}
				</li>
			{/each}
		</ul>
	</div>
</div>
<div class="row charts" use:dndzone={{ items: dsOrder }} onconsider={(evt) => (dsOrder = evt.detail.items)} onfinalize={(evt) => (dsOrder = evt.detail.items)}>
	{#each dsOrder as { id: dataset } (dataset)}
		<ChartDisplay {dataset} {config} {bigView} {bigViewCharts} />
	{/each}
</div>

<dialog
	open={!!(bigView && bigViewCharts[0] && bigViewCharts[1])}
	onclick={self(() => {
		bigView = false;
		bigViewCharts = [null, null];
	})}
>
	<article class="big-modal">
		{#if bigView && bigViewCharts[0] && bigViewCharts[1]}
			<div class="row">
				<div class="col center">
					{#await bigViewCharts[0]}
						Loading...
					{:then charts}
						<img src={charts.clustering} alt="Enlarged" />
						<img src={charts.violin} alt="Enlarged" />
					{/await}
				</div>
				<div class="col center">
					{#await bigViewCharts[1]}
						Loading...
					{:then charts}
						<img src={charts.clustering} alt="Enlarged" />
						<img src={charts.violin} alt="Enlarged" />
					{/await}
				</div>
			</div>
		{/if}
	</article>
</dialog>

<style lang="scss">
	@keyframes spin {
		0% {
			transform: rotate(0deg);
		}

		100% {
			transform: rotate(360deg);
		}
	}

	:global(.fa-spinner) {
		animation: spin 1s linear infinite;
	}

	.controls {
		height: 300px;
		gap: 1em;

		h2 {
			margin-bottom: 8px;
		}

		.gene-list {
			max-height: 180px;
			padding-block-start: 1px;
			padding-inline-start: 1px;
			padding-inline-end: 1px;
			overflow-y: scroll;

			.gene-option {
				list-style: none;
				padding: 0.2em 0.5em;
				cursor: pointer;
				opacity: 0.7;

				&:hover {
					opacity: 1;
				}

				&.selected {
					outline: 1px solid var(--primary-background);
				}
			}
		}
	}

	.charts {
		max-width: 100vw;
	}

	.big-modal {
		padding: 0;
		max-width: fit-content;

		.row {
			gap: 5vh;

			.col {
				gap: 5vh;
				width: 40vw;

				img {
					height: 40vh;
				}
			}
		}
	}
</style>
