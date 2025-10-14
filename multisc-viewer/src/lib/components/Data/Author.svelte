<script lang="ts">
	import type { Author } from '$lib/types/daemon';
	import { Portal, Tooltip } from '@skeletonlabs/skeleton-svelte';
	import Chip from '../List/Chip.svelte';

	let {
		author = [],
		tags = $bindable(),
		tag = 'author'
	}: { author?: Author; tags?: string[]; tag?: string } = $props();

	let open = $state(false);
</script>

<!-- onclick={() => {
		if (tags && tag && !tags.includes(tag)) {
			tags.push(`${tag}=${author[0]?.name}`);
		}
	}} -->
<Tooltip
	{open}
	onOpenChange={(e) => (open = e.open)}
	positioning={{ placement: 'left' }}
	openDelay={200}
	interactive
>
	<Tooltip.Trigger>
		<span class="chip preset-tonal">
			{author[0]?.name}{author.length > 1 ? ' et al.' : ''}
		</span>
	</Tooltip.Trigger>
	<Portal>
		<Tooltip.Positioner>
			<Tooltip.Content class="max-w-md card bg-surface-100-900 p-2 shadow-xl">
				<div class="flex flex-wrap gap-1">
					{#each author as a (`chip-a-${a.name}`)}
						<Chip bind:tags tag={`${tag}:${a.name}`}>
							{a.name}
						</Chip>
					{/each}
				</div>
			</Tooltip.Content>
		</Tooltip.Positioner>
	</Portal>
</Tooltip>
