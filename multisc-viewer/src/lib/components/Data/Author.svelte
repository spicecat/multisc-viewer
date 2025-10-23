<script lang="ts">
	import type { Author } from '$lib/types/daemon';
	import { Portal, Tooltip } from '@skeletonlabs/skeleton-svelte';
	import Chip from '../List/Chip.svelte';

	let {
		author = [],
		tags = $bindable(),
		tag = 'author'
	}: { author?: Author; tags?: string[]; tag?: string } = $props();
</script>

<Tooltip positioning={{ placement: 'left' }} openDelay={200} interactive>
	<Tooltip.Trigger
		type="button"
		onclick={(e) => {
			const tag = `author=${author[0]?.name}`;
			if (tags && !tags.includes(tag)) tags.push(tag);
			e.stopPropagation();
		}}
	>
		<span class="chip preset-tonal">{author[0]?.name}{author.length > 1 ? ', et al.' : ''}</span>
	</Tooltip.Trigger>
	<Portal>
		<Tooltip.Positioner>
			<Tooltip.Content class="max-w-md card bg-primary-100 p-2">
				<div class="inline-flex flex-wrap gap-1">
					{#each author as a (`chip-a-${a.name}`)}
						<Chip bind:tags tag={`${tag}=${a.name}`}>
							{a.name}
						</Chip>
					{/each}
				</div>
			</Tooltip.Content>
		</Tooltip.Positioner>
	</Portal>
</Tooltip>
