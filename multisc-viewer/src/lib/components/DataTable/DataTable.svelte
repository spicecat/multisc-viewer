<script lang="ts">
	import type { Columns, Data } from '$lib/types/data-table';
	import { ArrowRight } from '@lucide/svelte';
	import { Pagination, Progress } from '@skeletonlabs/skeleton-svelte';

	let {
		name,
		data,
		columns,
		select
	}: {
		name: string;
		data: Data[];
		columns: Columns;
		select?: string | string[];
	} = $props();

	let selected = $state<string | string[] | undefined>(select);

	let page = $state(1);
	let pageSize = $state(10);
	const slice = $derived(data.slice((page - 1) * pageSize, page * pageSize));

	$effect(() => {
		const totalPages = Math.ceil(data.length / pageSize);
		if (page > totalPages) page = Math.max(1, totalPages);
	});
</script>

<div class="table-wrap">
	<table class="table caption-bottom">
		<thead>
			<tr>
				{#if select}
					<th>
						<label class="flex items-center space-x-2">
							{#if Array.isArray(selected)}
								<input
									type="checkbox"
									class="checkbox"
									checked={selected.length === data.length}
									onchange={(e) => {
										selected = e.currentTarget.checked ? data.map((item) => item.id) : [];
									}}
								/>
								<p class="font-bold">{selected.length}</p>
							{:else}
								<p class="font-bold">{selected}</p>
							{/if}
						</label>
					</th>
				{/if}
				{#each columns as { label } (`label-${label}`)}
					<th>
						{label}
					</th>
				{/each}
			</tr>
		</thead>
		<tbody>
			{#await slice}
				<Progress value={null} />
			{:then}
				{#each slice as item (`item-${item.id}`)}
					{@const id = item.id}
					<tr>
						{#if select}
							<td>
								{#if Array.isArray(select)}
									<input type="checkbox" class="checkbox" {name} value={id} bind:group={selected} />
								{:else}
									<input type="radio" class="radio" {name} value={id} bind:group={selected} />
								{/if}
							</td>
						{/if}
						{#each columns as { key, href } (`column-${id}-${key}`)}
							{@const d = item[key]}
							<td>
								{#if href}
									<a class="btn preset-tonal-primary text-wrap" href={String(item[href])}>
										<span>{d}</span>
										<ArrowRight />
									</a>
								{:else if Array.isArray(d)}
									<div class="flex flex-wrap gap-1">
										{#each d as chip (`chip-${id}-${key}-${chip}`)}
											<span class="chip preset-tonal">{chip}</span>
										{/each}
									</div>
								{:else}
									<p>{d}</p>
								{/if}
							</td>
						{/each}
					</tr>
				{/each}
			{/await}
		</tbody>
	</table>
</div>
<footer class="flex justify-center">
	<Pagination {data} {page} onPageChange={(e) => (page = e.page)} {pageSize} alternative />
</footer>
