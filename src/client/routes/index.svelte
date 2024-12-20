<script lang="ts">
	import { onMount } from 'svelte';

	let open: boolean = false,
		datasets: string[] | null = null,
		selectedDatasets: string[] = [],
		genes: Promise<string[]> = Promise.resolve([]),
		geneSearch: string = '',
		selectedGene: string | null = null,
		groupBy: string = 'Genotype',
		plots: Record<string, { clustering: string; violin: string }> | null = null,
		ordering: string[] | null = null,
		plotting: boolean = false,
		addingDataset: boolean = false,
		newDatasets: FileList = null,
		newDatasetName: string = '',
		uploading: boolean = false;

	$: splitBy = groupBy === 'Genotype' ? 'CellType' : 'Genotype';

	let debounceTimeout: any | null = null;
	function fetchGenes(): void {
		if (debounceTimeout !== null) {
			clearTimeout(debounceTimeout);
		}

		debounceTimeout = setTimeout(() => {
			debounceTimeout = null;
			genes = fetch(`/genes?datasets=${selectedDatasets.join(',')}`).then((res) => res.json());
		}, 500);
		selectedGene = null;
	}

	function plot(): void {
		if (selectedGene !== null) {
			fetch(`/plots?datasets=${selectedDatasets.join(',')}&gene=${selectedGene}&groupBy=${groupBy}&splitBy=${splitBy}`)
				.then((res) => res.json())
				.then((plotsObj) => {
					plots = plotsObj;
					ordering = selectedDatasets;
					open = false;
					plotting = false;
				});
			plots = null;
			ordering = null;
			plotting = true;
		}
	}

	function upload(): void {
		if (newDatasets[0] && newDatasetName !== null) {
			const body = new FormData();

			body.append('data', newDatasets[0]);
			body.append('name', newDatasetName);

			console.log(newDatasets[0]);

			// fetch('/datasets', { body, method: 'POST' }).then(() => {
			// 	if (datasets) datasets = [...datasets, newDatasetName];
			// 	newDatasetName = '';
			// 	addingDataset = false;
			// 	uploading = false;
			// });
			uploading = true;
		}
	}

	onMount(() => {
		fetch('/datasets')
			.then((res) => res.json())
			.then((data) => (datasets = data));
	});
</script>

<main class="container">
	<h1>Dataset Comparison</h1>

	<button on:click={() => (open = true)}>Config</button>
	<button on:click={() => (addingDataset = true)}>Upload Dataset</button>

	{#if plots !== null && ordering !== null}
		<div class="row">
			{#each ordering as dataset}
				<div class="col">
					<h3 class="dataset">{dataset.replaceAll('_', ' ')}</h3>
					<img src={plots[dataset].clustering} alt="{dataset} clustering" />
					<img src={plots[dataset].violin} alt="{dataset} violin" />
				</div>
			{/each}
		</div>
	{/if}
</main>

<dialog {open} on:click|self={() => (open = false)}>
	<article class="modal">
		<header>
			<a href="/" class="close" on:click|preventDefault={() => (open = false)}></a>
			<h2>Config</h2>
		</header>
		<fieldset on:change={fetchGenes}>
			<legend>
				<h3>Datasets</h3>
			</legend>
			{#if datasets !== null}
				{#each datasets as dataset}
					<label>
						<input type="checkbox" value={dataset} bind:group={selectedDatasets} />
						{dataset}
					</label>
				{/each}
			{/if}
		</fieldset>
		<hr />
		<div class="row">
			<div class="col">
				<h4>Group By</h4>
				<p>{groupBy}</p>
			</div>
			<div class="col swap">
				<button on:click={() => (groupBy = splitBy)}>
					<i class="fa-solid fa-left-right"></i>
				</button>
			</div>
			<div class="col">
				<h4>Split By</h4>
				<p>{splitBy}</p>
			</div>
		</div>
		<hr />
		{#await genes}
			Loading genes...
		{:then genes}
			<label>
				Gene
				<input type="text" placeholder="Search" bind:value={geneSearch} />
			</label>
			{#if genes}
				{@const filteredGenes = genes.filter((gene) => gene.toLowerCase().includes(geneSearch.toLowerCase())).slice(0, 100)}
				<ul class="gene-list">
					{#each filteredGenes as gene}
						<li class="gene-option" class:selected={selectedGene === gene} on:click={() => (selectedGene = gene)}>{gene}</li>
					{/each}
				</ul>
				{#if genes.length > filteredGenes.length}
					<p>{genes.length - filteredGenes.length} more...</p>
				{/if}
			{/if}
		{:catch err}
			Error: {err}
		{/await}
		<footer>
			<button on:click={plot} disabled={plotting}>
				{#if plotting}
					<i class="fa-solid fa-spinner"></i>
				{/if}
				Plot
			</button>
		</footer>
	</article>
</dialog>
<dialog open={addingDataset} on:click|self={() => (addingDataset = false)}>
	<article>
		<header>
			<a href="/" class="close" on:click|preventDefault={() => (addingDataset = false)}></a>
			<h2>Upload Dataset</h2>
		</header>
		<h3>RDS File</h3>
		<input type="file" bind:files={newDatasets} />
		<label>
			Name
			<input type="text" bind:value={newDatasetName} />
		</label>
		<footer>
			<button on:click={upload} disabled={uploading}>
				{#if uploading}
					<i class="fa-solid fa-spinner"></i>
				{/if}
				Upload
			</button>
		</footer>
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

	img {
		width: 300px;
	}

	.dataset {
		text-align: center;
	}

	.swap {
		justify-content: flex-end;
	}

	.modal {
		height: 800px;
		overflow-y: scroll;

		.gene-list {
			max-height: 400px;
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

	:global(.fa-spinner) {
		animation: spin 1s linear infinite;
	}
</style>
