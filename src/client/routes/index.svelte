<script lang="ts">
	// import { goto } from "$app/navigation";
	import DataTable from '$lib/components/DataTable.svelte';
	import Button, { Label } from '@smui/button';
	import IconButton from '@smui/icon-button';
	import TextField from '@smui/textfield';
	import { onMount } from 'svelte';

	const { token }: IndexProps = $props();

	let datasets: object[] = $state([]),
		selectedDatasets: string[] = $state([]),
		datasetSearch: string = $state(''),
		filteredDatasets: object[] = $derived(filterItems(datasets, datasetSearch)),
		genes: object[] = $state([]),
		selectedGene: string = $state(''),
		geneSearch: string = $state(''),
		filteredGenes: object[] = $derived(filterItems(genes, geneSearch)),
		groupBy: string = $state('Genotype'),
		splitBy = $derived(groupBy === 'Genotype' ? 'CellType' : 'Genotype');

	// TODO: add lodash debounce
	$effect(() => {
		if (selectedDatasets.length) {
			fetch(`/genes?datasets=${selectedDatasets.join(',')}`)
				.then((res) => res.json())
				.then((data) => (genes = data.map((gene: string) => ({ gene }))));
			fetch(`/preload?token=${token}&datasets=${selectedDatasets.join(',')}`, { method: 'POST' });
		} else {
			genes = [];
			selectedGene = '';
		}
	});

	function plot(): void {
		if (selectedGene && selectedDatasets.length)
			window.location.href = `/compare?datasets=${selectedDatasets.join(',')}&gene=${selectedGene}&groupBy=${groupBy}&splitBy=${splitBy}`; // use goto
	}

	function filterItems(items: object[], search: string): object[] {
		return items.filter((item) => Object.values(item).some((value) => value.toString().toLowerCase().includes(search.toLowerCase())));
	}

	onMount(() => {
		fetch('/datasets')
			.then((res) => res.json())
			.then((data) => (datasets = data));
	});

	const datasetColumns = [
			{ key: 'name', label: 'Name' },
			{ key: 'year', label: 'Year' },
			{ key: 'region', label: 'Region' },
			{ key: 'PMID', label: 'PMID' },
			{ key: 'species', label: 'Species' },
			{ key: 'author', label: 'Author' },
			{ key: 'disease', label: 'Disease' }
		],
		geneColumns = [{ key: 'gene', label: 'Gene' }];
</script>

<svelte:head>
	<!-- SMUI -->
	<link rel="stylesheet" href="https://unpkg.com/svelte-material-ui/themes/unity.css" />
</svelte:head>

<main class="container">
	<h1>Dataset Comparison</h1>

	<div class="row">
		<div>
			<h4>Group By:</h4>
			<p>{groupBy}</p>
		</div>
		<IconButton
			class="material-icons"
			onclick={() => {
				groupBy = splitBy;
			}}>swap_horiz</IconButton
		>
		<div>
			<h4>Split By:</h4>
			<p>{splitBy}</p>
		</div>
		<Button onclick={plot} variant="raised">
			<Label>Plot</Label>
		</Button>
	</div>

	<TextField bind:value={datasetSearch} label="Datasets">
		{#snippet leadingIcon()}
			<IconButton class="material-icons">search</IconButton>
		{/snippet}
	</TextField>
	<DataTable columns={datasetColumns} data={filteredDatasets} bind:selected={selectedDatasets} />
	<TextField bind:value={geneSearch} label="Genes">
		{#snippet leadingIcon()}
			<IconButton class="material-icons">search</IconButton>
		{/snippet}
	</TextField>
	<DataTable columns={geneColumns} data={filteredGenes} bind:selected={selectedGene} />
</main>

<style>
	.row {
		display: flex;
		align-items: center;
	}
	.row > div {
		margin-right: 1em;
	}
</style>
