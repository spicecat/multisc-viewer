<script lang="ts">
	import type { Columns, Data, Select } from '$lib/types/data-table';
	import { Pagination } from '@skeletonlabs/skeleton-svelte';

	let {
		data,
		columns,
		select
	}: {
		data: Data[];
		columns: Columns;
		select: Select;
	} = $props();

	let selected = $state<string | string[]>(select === 'checkbox' ? [] : '');

	let page = $state(1);
	let pageSize = $state(10);
	const slice = $derived(data.slice((page - 1) * pageSize, page * pageSize));
</script>

<section>
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
									name="selected"
									onchange={(e) => {
										selected = e.currentTarget.checked
											? data.map((item) => (typeof item === 'object' ? item.id : item))
											: [];
									}}
								/>
							{/if}
						</th>
					{/if}
					{#each columns as { label }}
						<th>
							{label}
						</th>
					{/each}
				</tr>
			</thead>
			<tbody>
				{#each slice as item}
					{@const id = typeof item === 'object' ? item.id : item}
					<tr>
						{#if select}
							<td>
								{#if select === 'checkbox'}
									<input
										type="checkbox"
										class="checkbox"
										name="selected"
										value={id}
										bind:group={selected}
									/>
								{:else}
									<input
										type="radio"
										class="radio"
										name="selected"
										value={id}
										bind:group={selected}
									/>
								{/if}
							</td>
						{/if}
						{#each columns as { key, href }}
							{@const d =
								typeof key === 'function' ? key(item) : typeof item === 'object' ? item[key] : item}
							<td>
								{#if href}
									<a class="btn preset-tonal text-wrap" href={href(item)}>{d}</a>
								{:else if Array.isArray(d)}
									<div class="flex flex-wrap gap-1">
										{#each d as chip}
											<span class="chip preset-tonal">{chip}</span>
										{/each}
									</div>
								{:else}
									{d}
								{/if}
							</td>
						{/each}
					</tr>
				{/each}
			</tbody>
		</table>
	</div>
	<footer class="flex justify-center">
		<Pagination {data} {page} onPageChange={(e) => (page = e.page)} {pageSize} alternative />
	</footer>
</section>
