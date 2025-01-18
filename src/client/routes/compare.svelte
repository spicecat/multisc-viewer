<script lang="ts">
	import meta from '$meta';

	const { genes, gene }: CompareProps = $props();

	const { query } = meta;
	const ordering = query.datasets.split(',');

	let geneSearch: string = $state(''),
		groupBy: string = $state('Genotype'),
		selectedGene: string = $state(gene);

	const filteredGenes = $derived(genes.filter((gene) => gene.toLowerCase().includes(geneSearch.toLowerCase())).slice(0, 100)),
		splitBy = $derived(groupBy === 'Genotype' ? 'CellType' : 'Genotype');

	$effect(() => {
		if (!history.state || history.state.selectedGene !== selectedGene || history.state.groupBy !== groupBy) {
			history.pushState({ selectedGene, groupBy }, '', `/compare?datasets=${query.datasets}&gene=${selectedGene}&groupBy=${groupBy}&splitBy=${splitBy}`);
		}
	});
</script>

<div class="row controls">
	<div class="gene-select">
		<label>
			<h2>Gene</h2>
			<input type="text" placeholder="Search" bind:value={geneSearch} />
		</label>
		<ul class="gene-list">
			{#each filteredGenes as gene}
				<li class="gene-option" class:selected={selectedGene === gene} onclick={() => (selectedGene = gene)}>{gene}</li>
			{/each}
		</ul>
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
</div>
{#await fetch(`/plots?datasets=${query.datasets}&gene=${selectedGene}&groupBy=${groupBy}&splitBy=${splitBy}`).then((res) => res.json())}
	<i class="fa-solid fa-spinner"></i> Generating Plots...
{:then plots}
	<div class="row">
		{#each ordering as dataset}
			<div class="col">
				<h3 class="dataset">{dataset.replaceAll('_', ' ')}</h3>
				<img src={plots[dataset].clustering} alt="{dataset} clustering" />
				<img src={plots[dataset].violin} alt="{dataset} violin" />
			</div>
		{/each}
	</div>
{:catch err}
	<h1>Error: {err}</h1>
{/await}

<style lang="scss">
	@keyframes spin {
		0% {
			transform: rotate(0deg);
		}

		100% {
			transform: rotate(360deg);
		}
	}

	img {
		width: 300px;
	}

	.dataset {
		text-align: center;
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
</style>
