<script lang="ts">
	import { page } from '$app/state';
	import type { Datasets, DEGs, Genes } from '$lib/types/daemon';
	import { Dna } from '@lucide/svelte';
	import { SvelteMap } from 'svelte/reactivity';
	import DataTable from '../List/DataTable.svelte';
	import Search from '../List/Search.svelte';

	let { datasets, genes, degs }: { datasets: Datasets; genes: Genes; degs: DEGs } = $props();
	let data = $derived.by(() => {
		const genesMap = new SvelteMap<string, { _id: string; genes: string[]; degs: string[] }>();
		for (const [ds, dsGenes] of Object.entries(genes)) {
			for (const gene of dsGenes) {
				if (!genesMap.has(gene)) genesMap.set(gene, { _id: gene, genes: [], degs: [] });
				genesMap.get(gene)?.genes.push(ds);
			}
		}
		for (const [ds, dsDEGs] of Object.entries(degs)) {
			for (const gene of dsDEGs) {
				if (!genesMap.has(gene)) genesMap.set(gene, { _id: gene, genes: [], degs: [] });
				genesMap.get(gene)?.degs.push(ds);
			}
		}
		return Array.from(genesMap.values()).sort(
			(a, b) =>
				b.degs.length - a.degs.length ||
				b.genes.length - a.genes.length ||
				a._id.localeCompare(b._id)
		);
	});
	let items: typeof data = $state([]);

	const searchOptions = {
		keys: ['_id', 'genes', 'degs']
	};

	let tags: string[] = $state([]);
	let selected = $state(page.url.searchParams.getAll('gene'));
</script>

<Search name="Genes" {data} bind:items bind:tags {searchOptions} />
<DataTable data={items} selected={data.filter((g) => selected.includes(g._id))}>
	<tr>
		<th>
			<label class="flex items-center space-x-2">
				<input
					type="checkbox"
					class="checkbox"
					checked={items.every((gene) => selected.includes(gene._id))}
					onchange={(e) =>
						(selected = e.currentTarget.checked ? items.map((gene) => gene._id) : [])}
				/>
				<p class="font-bold">{selected.length}</p>
			</label>
		</th>
		<th>Gene</th>
		{#each Object.keys(genes) as ds (`genes-ds-${ds}`)}
			{@const name = datasets[ds]?.displayName ?? datasets[ds]?.name ?? ds}
			<th>
				{name}
				<p class="font-bold">
					{genes[ds]?.length} genes,
					<span class="text-red-500">{degs[ds]?.length} DEGs</span>
				</p>
			</th>
		{/each}
		<th>DEGs</th>
	</tr>
	{#snippet row(gene)}
		<td>
			<input type="checkbox" class="checkbox" name="gene" value={gene._id} bind:group={selected} />
		</td>
		<td>{gene._id}</td>
		{#each Object.keys(genes) as ds (`genes-ds-${ds}`)}
			<td>
				{#if gene.degs.includes(ds)}
					<Dna color="red" />
				{:else if gene.genes.includes(ds)}
					<Dna />
				{/if}
			</td>
		{/each}
		<td>{gene.degs.length}</td>
	{/snippet}
</DataTable>

<!-- Prevent resetting selected on page change -->
{#each selected as s, i (`gene-hidden-${i}`)}
	<input type="checkbox" class="hidden" value={s} bind:group={selected} />
{/each}
