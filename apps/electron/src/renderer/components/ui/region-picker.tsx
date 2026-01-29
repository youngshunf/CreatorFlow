/**
 * RegionPicker - Cascading region selection component.
 *
 * Three-level cascading select for province, city, and district.
 */

import * as React from 'react'
import { ChevronDown } from 'lucide-react'

import { cn } from '@/lib/utils'
import { Button } from '@/components/ui/button'
import {
  Popover,
  PopoverContent,
  PopoverTrigger,
} from '@/components/ui/popover'
import { ScrollArea } from '@/components/ui/scroll-area'
import {
  getProvinces,
  getCitiesByProvince,
  getDistrictsByCity,
  type Region,
} from '@/data/china-regions'

export interface RegionValue {
  province?: string
  city?: string
  district?: string
}

export interface RegionPickerProps {
  /** Current value */
  value?: RegionValue
  /** Callback when region changes */
  onChange?: (value: RegionValue) => void
  /** Placeholder text */
  placeholder?: string
  /** Disabled state */
  disabled?: boolean
  /** Additional className */
  className?: string
}

type Tab = 'province' | 'city' | 'district'

/**
 * RegionPicker - Cascading province/city/district selector
 *
 * @example
 * <RegionPicker
 *   value={{ province: '广东省', city: '深圳市', district: '南山区' }}
 *   onChange={setRegion}
 * />
 */
export function RegionPicker({
  value,
  onChange,
  placeholder = '请选择地区',
  disabled = false,
  className,
}: RegionPickerProps) {
  const [open, setOpen] = React.useState(false)
  const [activeTab, setActiveTab] = React.useState<Tab>('province')

  // Get available options based on current selection
  const provinces = React.useMemo(() => getProvinces(), [])
  const cities = React.useMemo(
    () => (value?.province ? getCitiesByProvince(value.province) : []),
    [value?.province]
  )
  const districts = React.useMemo(
    () =>
      value?.province && value?.city
        ? getDistrictsByCity(value.province, value.city)
        : [],
    [value?.province, value?.city]
  )

  // Display value
  const displayValue = React.useMemo(() => {
    const parts = [value?.province, value?.city, value?.district].filter(Boolean)
    return parts.length > 0 ? parts.join(' / ') : null
  }, [value])

  // Handle province selection
  const handleProvinceSelect = React.useCallback(
    (region: Region) => {
      onChange?.({ province: region.name, city: undefined, district: undefined })
      setActiveTab('city')
    },
    [onChange]
  )

  // Handle city selection
  const handleCitySelect = React.useCallback(
    (region: Region) => {
      onChange?.({ ...value, city: region.name, district: undefined })
      setActiveTab('district')
    },
    [onChange, value]
  )

  // Handle district selection
  const handleDistrictSelect = React.useCallback(
    (region: Region) => {
      onChange?.({ ...value, district: region.name })
      setOpen(false)
    },
    [onChange, value]
  )

  // Reset to province tab when opening
  React.useEffect(() => {
    if (open) {
      if (value?.city && value?.province) {
        setActiveTab('district')
      } else if (value?.province) {
        setActiveTab('city')
      } else {
        setActiveTab('province')
      }
    }
  }, [open, value])

  return (
    <Popover open={open} onOpenChange={setOpen}>
      <PopoverTrigger asChild>
        <Button
          variant="outline"
          disabled={disabled}
          className={cn(
            'w-full justify-between text-left font-normal',
            !displayValue && 'text-muted-foreground',
            className
          )}
        >
          <span className="truncate">{displayValue || placeholder}</span>
          <ChevronDown className="ml-2 size-4 shrink-0 opacity-50" />
        </Button>
      </PopoverTrigger>
      <PopoverContent className="w-[320px] p-0" align="start">
        {/* Tabs */}
        <div className="flex border-b">
          <TabButton
            active={activeTab === 'province'}
            onClick={() => setActiveTab('province')}
            disabled={false}
          >
            {value?.province || '省份'}
          </TabButton>
          <TabButton
            active={activeTab === 'city'}
            onClick={() => setActiveTab('city')}
            disabled={!value?.province}
          >
            {value?.city || '城市'}
          </TabButton>
          <TabButton
            active={activeTab === 'district'}
            onClick={() => setActiveTab('district')}
            disabled={!value?.city}
          >
            {value?.district || '区县'}
          </TabButton>
        </div>

        {/* Content */}
        <ScrollArea className="h-[280px]">
          <div className="p-2">
            {activeTab === 'province' && (
              <RegionGrid
                regions={provinces}
                selectedName={value?.province}
                onSelect={handleProvinceSelect}
              />
            )}
            {activeTab === 'city' && (
              <RegionGrid
                regions={cities}
                selectedName={value?.city}
                onSelect={handleCitySelect}
              />
            )}
            {activeTab === 'district' && (
              <RegionGrid
                regions={districts}
                selectedName={value?.district}
                onSelect={handleDistrictSelect}
              />
            )}
          </div>
        </ScrollArea>
      </PopoverContent>
    </Popover>
  )
}

interface TabButtonProps {
  children: React.ReactNode
  active: boolean
  disabled: boolean
  onClick: () => void
}

function TabButton({ children, active, disabled, onClick }: TabButtonProps) {
  return (
    <button
      type="button"
      onClick={onClick}
      disabled={disabled}
      className={cn(
        'flex-1 px-3 py-2 text-sm font-medium transition-colors truncate',
        'border-b-2 -mb-px',
        active
          ? 'border-primary text-primary'
          : 'border-transparent text-muted-foreground hover:text-foreground',
        disabled && 'opacity-50 cursor-not-allowed'
      )}
    >
      {children}
    </button>
  )
}

interface RegionGridProps {
  regions: Region[]
  selectedName?: string
  onSelect: (region: Region) => void
}

function RegionGrid({ regions, selectedName, onSelect }: RegionGridProps) {
  if (regions.length === 0) {
    return (
      <div className="py-8 text-center text-sm text-muted-foreground">
        请先选择上一级
      </div>
    )
  }

  return (
    <div className="grid grid-cols-3 gap-1">
      {regions.map((region) => (
        <button
          key={region.code}
          type="button"
          onClick={() => onSelect(region)}
          className={cn(
            'px-2 py-1.5 text-sm rounded-md truncate text-left transition-colors',
            'hover:bg-muted',
            selectedName === region.name && 'bg-primary/10 text-primary font-medium'
          )}
        >
          {region.name}
        </button>
      ))}
    </div>
  )
}
