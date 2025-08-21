<script lang="ts">
	import type { Columns, Data, Select } from '$lib/types/data-table';
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
		select?: Select;
	} = $props();

	let selected = $state<string | string[]>(select === 'checkbox' ? [] : '');

	let page = $state(1);
	let pageSize = $state(10);
	const slice = $derived(data.slice((page - 1) * pageSize, page * pageSize));
</script>

<div class="table-wrap">
	<table class="table caption-bottom">
		<thead>
			<tr>
				{#if select}
					<th>
						{#if select === 'checkbox'}
							<input
								type="checkbox"
								class="checkbox"
								checked={selected.length === data.length}
								onchange={(e) => {
									selected = e.currentTarget.checked ? data.map((item) => item.id) : [];
								}}
							/>
						{/if}
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
					{@const id = typeof item === 'object' ? item.id : item}
					<tr>
						{#if select}
							<td>
								{#if select === 'checkbox'}
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
									<a class="btn preset-tonal text-wrap" href={String(item.href)}>{d} &rarr;</a>
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
