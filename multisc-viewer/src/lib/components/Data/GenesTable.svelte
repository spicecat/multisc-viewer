<script lang="ts">
	import { page } from '$app/state';
	import type { Datasets, GenesRows } from '$lib/types/daemon';
	import { Dna } from '@lucide/svelte';
	import DataTable from '../List/DataTable.svelte';
	import Search from '../List/Search.svelte';

	let { datasets, genesRows }: { datasets: Datasets; genesRows: GenesRows } = $props();

	let selected = $state(page.url.searchParams.getAll('gene'));

	let items: GenesRows = $state([]);
	let tags: string[] = $state([]);

	const searchOptions = { keys: ['_id', 'datasets', 'degs'] };
</script>

<Search
	name="Differentially Expressed Genes"
	data={genesRows}
	bind:items
	bind:tags
	{searchOptions}
/>
<DataTable data={items} selected={genesRows.filter((g) => selected.includes(g._id))}>
	<tr>
		<th>
			<label class="flex items-center space-x-2">
				<input
					type="checkbox"
					class="checkbox"
					aria-label="Select all genes"
					checked={items.every((gene) => selected.includes(gene._id))}
					onchange={(e) =>
						(selected = e.currentTarget.checked ? items.map((gene) => gene._id) : [])}
				/>
				<p class="font-bold">{selected.length}</p>
			</label>
		</th>
		<th>Gene</th>
		{#each Object.keys(datasets) as ds (`genes-th-ds-${ds}`)}
			<th class="border-x border-gray-300">
				{datasets[ds]?.displayName ?? ds}
			</th>
		{/each}
		<th># of Datasets</th>
	</tr>
	{#snippet row(gene)}
		<td>
			<input
				type="checkbox"
				class="checkbox"
				aria-label="Select gene {gene._id}"
				name="gene"
				value={gene._id}
				bind:group={selected}
			/>
		</td>
		<td>{gene._id}</td>
		{#each Object.keys(datasets) as ds (`genes-td-ds-${ds}`)}
			<td class="border-x border-gray-300">
				<div class="flex justify-center">
					{#if gene.datasets.includes(ds)}
						<Dna />
					{/if}
				</div>
			</td>
		{/each}
		<td class="text-center">{gene.datasets.length}</td>
	{/snippet}
</DataTable>

<!-- Prevent resetting selected on page change -->
{#each selected as s, i (`gene-hidden-${i}`)}
	<input type="checkbox" class="hidden" value={s} bind:group={selected} />
{/each}
