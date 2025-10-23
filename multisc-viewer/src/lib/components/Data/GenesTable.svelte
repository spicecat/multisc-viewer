<script lang="ts">
	import { page } from '$app/state';
	import type { Datasets, DEGs, Genes } from '$lib/types/daemon';
	import { Dna } from '@lucide/svelte';
	import { SvelteMap, SvelteSet } from 'svelte/reactivity';
	import DataTable from '../List/DataTable.svelte';
	import Search from '../List/Search.svelte';

	let { datasets, genes, degs }: { datasets: Datasets; genes: Genes; degs: DEGs } = $props();

	let selected = $state(page.url.searchParams.getAll('gene'));

	let data = $derived.by(() => {
		const genesMap = new SvelteMap<
			string,
			{ _id: string; datasets: SvelteSet<string>; degs: SvelteSet<string> }
		>();
		for (const [ds, dsDEGs] of Object.entries(degs))
			for (const [deg, dsDEG] of Object.entries(dsDEGs))
				for (const g of dsDEG.gene) {
					if (!genesMap.has(g))
						genesMap.set(g, { _id: g, datasets: new SvelteSet(), degs: new SvelteSet() });
					genesMap.get(g)!.degs.add(deg);
					genesMap.get(g)!.datasets.add(ds);
				}

		return Array.from(genesMap.values()).sort(
			(a, b) => b.degs.size - a.degs.size || a._id.localeCompare(b._id)
		);
	});
	let items: typeof data = $state([]);
	let tags: string[] = $state([]);

	const searchOptions = { keys: ['_id'] };
</script>

<Search name="Differentially Expressed Genes" {data} bind:items bind:tags {searchOptions} />
<DataTable data={items} selected={data.filter((g) => selected.includes(g._id))}>
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
		{#each Object.keys(genes) as ds (`genes-th-ds-${ds}`)}
			{#if degs[ds]}
				{#each Object.keys(degs[ds]) as deg, i (`degs-th-deg-${deg}`)}
					{@const DEG = degs[ds][deg]}
					<th class={i === 0 ? 'border-l border-gray-300' : undefined}>
						{DEG.name ?? DEG._id}
						<p class="font-bold">
							<span>{genes[ds]?.length} genes</span>
							<span class="text-red-500">{DEG.gene.length} DEGs</span>
						</p>
					</th>
				{/each}
			{:else}
				<th class="border-x border-gray-300">
					{datasets[ds]?.displayName ?? ds}
					<p class="font-bold">
						<span>{genes[ds]?.length} genes</span>
					</p>
				</th>
			{/if}
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
		{#each Object.keys(genes) as ds (`genes-td-ds-${ds}`)}
			{#if degs[ds]}
				{#each Object.keys(degs[ds]) as deg, i (`degs-td-deg-${deg}`)}
					<td class={i === 0 ? 'border-l border-gray-300' : undefined}>
						<div class="flex justify-center">
							{#if gene.degs.has(deg)}
								<Dna color="red" />
							{:else if genes[ds].includes(gene._id)}
								<Dna />
							{/if}
						</div>
					</td>
				{/each}
			{:else}
				<td class="border-x border-gray-300">
					<div class="flex justify-center">
						{#if genes[ds].includes(gene._id)}
							<Dna />
						{/if}
					</div>
				</td>
			{/if}
		{/each}
		<td class="text-center">{gene.datasets.size}</td>
	{/snippet}
</DataTable>

<!-- Prevent resetting selected on page change -->
{#each selected as s, i (`gene-hidden-${i}`)}
	<input type="checkbox" class="hidden" value={s} bind:group={selected} />
{/each}
